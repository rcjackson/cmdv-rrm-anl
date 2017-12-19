import matplotlib
matplotlib.use('Agg')
import pyart
from matplotlib import pyplot as plt
import numpy as np
import glob
import os
from copy import deepcopy
from distributed import Client

import dask.bag as db
from time import sleep
from datetime import timedelta, datetime
from scipy.optimize import fmin_l_bfgs_b
import time
import time_procedures
import sys

sys.path.append('/home/rjackson/cmdv-rrm-anl/code/')

# File paths
berr_data_file_path = '/lcrc/group/earthscience/radar/stage/radar_disk_two/berr_rapic/'
data_path_cpol = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_rapic/'
out_file_path = '/lcrc/group/earthscience/rjackson/quicklook_plots/berr/'

## Berrima - 2009-2011 (new format), 2005-2005 (old format)
start_year = int(sys.argv[1])
start_month = int(sys.argv[2])
start_day = int(sys.argv[3])
start_hour = 0
start_minute = 1 

end_year = int(sys.argv[4])
end_month = int(sys.argv[5])
end_day = int(sys.argv[6])
end_hour = 0
end_minute = 2

def display_time(rad_date):
    import sys
    sys.path.append('/home/rjackson/cmdv-rrm-anl/code/')
    
    # Get a Radar object given a time period in the CPOL dataset
    data_path_cpol = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_rapic/'
    out_file_path = '/lcrc/group/earthscience/rjackson/quicklook_plots/berr/'
    out_data_path = '/lcrc/group/earthscience/rjackson/berr/'

    # CPOL in lassen or rapic?
    cpol_format = 1    # 0 = lassen, 1 = rapic
    
    def get_radar_from_cpol_rapic(time):
        from datetime import timedelta, datetime
        year_str = "%04d" % time.year
        month_str = "%02d" % time.month
        day_str = "%02d" % time.day
        hour_str = "%02d" % time.hour
        minute_str = "%02d" % time.minute
        second_str = "%02d" % time.second
        if(time.year == 2009 or (time.year == 2010 and
                                 time.month < 6)):
            dir_str = 'cpol_0910/rapic/'
        else:
            dir_str = 'cpol_1011/rapic/'

        file_name_str = (data_path_cpol +
                         dir_str + 
                         year_str +
                         month_str +
                         day_str +
                         hour_str +
                         minute_str +
                        'Gunn_Pt' +
                        '.rapic')
        radar = pyart.aux_io.read_radx(file_name_str)
        return radar 
    
    one_day_later = rad_date+timedelta(days=1)
    times, dates = time_procedures.get_radar_times_berr_cfradial(rad_date.year, 
                                                                 rad_date.month,
                                                                 rad_date.day,
                                                                 0, 
                                                                 1,
                                                                 one_day_later.year, 
                                                                 one_day_later.month,
                                                                 one_day_later.day,
                                                                 0, 
                                                                 2,
                                                                 )
    print(times)
    for rad_time in times:
        year_str = "%04d" % rad_time.year
        month_str = "%02d" % rad_time.month
        day_str = "%02d" % rad_time.day
        hour_str = "%02d" % rad_time.hour
        minute_str = "%02d" % rad_time.minute
        second_str = "%02d" % rad_time.second

        # Check to see if Cf/Radial file already exists...
        out_path = (out_data_path +
	           '/' +
                    year_str +
                    '/' +
	            month_str +
                    '/' +
                    day_str +
	            '/')
        if not os.path.exists(out_path):
            os.makedirs(out_path)

        out_file = ('BerrimaVol' + 
                    year_str +
                    month_str +
                    day_str + 
                    hour_str + 
                    minute_str + 
                    second_str +
                    '_deal.cf')
    
        if(not os.path.isfile(out_file)):
            try:
                radar = time_procedures.get_radar_from_berr_cfradial(rad_time)
                if(not 'last_Radar' in locals()):
                    last_Radar = radar
                if(cpol_format == 1):
                    ref_field = 'Refl'
                    vel_field = 'Vel'
                else:
                    ref_field = 'reflectivity'
                    vel_field = 'velocity'

                # Get sounding for 4DD intialization
                one_day_ago = rad_time-timedelta(days=1, minutes=1)
                sounding_times = time_procedures.get_sounding_times(one_day_ago.year,
                                                                    one_day_ago.month,
                                                                    one_day_ago.day,
                                                                    one_day_ago.hour,
                                                                    one_day_ago.minute,
                                                                    rad_time.year,
                                                                    rad_time.month,
                                                                    rad_time.day,
                                                                    rad_time.hour,
                                                                    rad_time.minute,
                                                                    minute_interval=60)
             
                sounding_time = sounding_times[len(sounding_times)-1]
                Sounding_netcdf = time_procedures.get_sounding(sounding_time)

                # Convert timestamps to datetime format
                Time = Sounding_netcdf.variables['time_offset'][:]
                base_time = Sounding_netcdf.variables['base_time'][:]
                alt = Sounding_netcdf.variables['alt'][:]
                u = Sounding_netcdf.variables['u_wind'][:]
                v = Sounding_netcdf.variables['v_wind'][:]
    
                Sounding_netcdf.close()
                steps = int(np.floor(len(u)/70))
                wind_profile = pyart.core.HorizontalWindProfile.from_u_and_v(alt[0::steps],
                                                                             u[0::steps],
                                                                             v[0::steps])
                
                ## 4DD expects speed, direction but HorizontalWindProfile outputs u_wind, v_wind
                wind_profile.u = wind_profile.u_wind
                wind_profile.v = wind_profile.v_wind
                sim_vel = pyart.util.simulated_vel_from_profile(radar, wind_profile, sim_vel_field=vel_field)
                radar.add_field('sim_velocity', sim_vel, replace_existing=True)

                # Dealias velocities
                gatefilter = pyart.correct.GateFilter(radar)
                
                gatefilter.exclude_below(ref_field, 0)
                gatefilter.exclude_masked(vel_field)
                gatefilter.exclude_invalid(vel_field)
                gatefilter.exclude_masked(ref_field)
                gatefilter.exclude_invalid(ref_field)
                gatefilter.exclude_above(ref_field, 80)
                gatefilter.exclude_below(vel_field, -75)
                gatefilter.exclude_above(vel_field, 75)
                
                texture = pyart.retrieve.calculate_velocity_texture(radar, vel_field=vel_field, wind_size=4)
                radar.add_field('velocity_texture', texture, replace_existing=True)
                gatefilter.exclude_above('velocity_texture', 2)
                gatefilter = pyart.correct.despeckle.despeckle_field(radar,
                                                                     vel_field,
                                                                     gatefilter=gatefilter,
                                                                     size=25) 

                # Calculate result from region based dealiasing
                corrected_velocity_region = pyart.correct.dealias_region_based(radar,
                                                                               vel_field=vel_field,
                                                                               keep_original=False,
                                                                               centered=True,
                                                                               gatefitler=gatefilter,
                                                                               interval_splits=6,
                                                                               skip_between_rays=2000,
                                                                               skip_along_ray=2000,
                                                                               rays_wrap_around=True,
                                                                               valid_min=-75,
                                                                               valid_max=75)   
                gfilter = gatefilter.gate_excluded
                vels = deepcopy(corrected_velocity_region['data'])
                vels_uncorr = radar.fields[vel_field]['data']
                sim_vels = radar.fields['sim_velocity']['data']
                v_nyq_vel = radar.instrument_parameters['nyquist_velocity']['data'][0]
                region_means = []
                regions = np.zeros(vels.shape)
                for nsweep, sweep_slice in enumerate(radar.iter_slice()):
                    sfilter = gfilter[sweep_slice]
                    vels_slice = vels[sweep_slice]
                    svels_slice = sim_vels[sweep_slice]
                    vels_uncorrs = vels_uncorr[sweep_slice]
                    valid_sdata = vels_uncorrs[~sfilter]
                    int_splits = pyart.correct.region_dealias._find_sweep_interval_splits(v_nyq_vel, 3, valid_sdata, nsweep)
                    regions[sweep_slice], nfeatures = pyart.correct.region_dealias._find_regions(vels_uncorrs, sfilter, limits=int_splits)


	            ## Minimize cost function that is sum of difference between regions and
                    def cost_function(nyq_vector):
                        cost = 0
                        i = 0
                        for reg in np.unique(regions[sweep_slice]):
                            add_value = np.abs(np.ma.mean(vels_slice[regions[sweep_slice] == reg]) + nyq_vector[i]*2*v_nyq_vel 
                                - np.ma.mean(svels_slice[regions[sweep_slice] == reg])) 
                
                            if(np.isfinite(add_value)):
                                cost += add_value
                            i = i + 1
                        return cost

	    
                    def gradient(nyq_vector):
                        gradient_vector = np.zeros(len(nyq_vector))
                        i = 0
                        for reg in np.unique(regions[sweep_slice]):
                            add_value = (np.ma.mean(vels_slice[regions[sweep_slice] == reg]) + nyq_vector[i]*2*v_nyq_vel
                                - np.ma.mean(svels_slice[regions[sweep_slice] == reg])) 
                            if(add_value > 0):
                                gradient_vector[i] = 2*v_nyq_vel
                            else:
                                gradient_vector[i] = -2*v_nyq_vel
                            i = i + 1
                        return gradient_vector


                    bounds_list = [(x,y) for (x,y) in zip(-5*np.ones(nfeatures+1), 5*np.ones(nfeatures+1))]
                    nyq_adjustments = fmin_l_bfgs_b(cost_function, np.zeros((nfeatures+1)), disp=True, fprime=gradient,
                                                    bounds=bounds_list, maxiter=30)
                    i = 0
                    for reg in np.unique(regions[sweep_slice]):
                        vels_slice[regions[sweep_slice] == reg] += v_nyq_vel*np.round(nyq_adjustments[0][i])
                        i = i + 1
                        vels[sweep_slice] = vels_slice

                corrected_velocity_region['data'] = vels                                           
                      
                radar.add_field_like(vel_field, 
                                     'corrected_velocity', 
                                     corrected_velocity_region['data'],
                                     replace_existing=True)


                time_procedures.write_radar_to_berr_cfradial(radar, rad_time)
                last_Radar = radar
                out_path = (out_file_path +
                            '/' +
                            year_str +
                            '/' +
	                    month_str +
                            '/' +
                            day_str +
                            '/')
                if not os.path.exists(out_path):
                    os.makedirs(out_path)

                out_file = hour_str + minute_str + '.png'
                plt.figure(figsize=(7,21))
                plt.subplot(411)

                display = pyart.graph.RadarMapDisplay(radar)
                display.plot_ppi(ref_field, 
                                 sweep=0, 
                                 cmap=pyart.graph.cm.NWSRef,
	                         vmin=0, 
                                 vmax=70,
                                 gatefilter=gatefilter)
                plt.subplot(412)
                display.plot_ppi('corrected_velocity', 
                                 sweep=0,
	                         cmap=pyart.graph.cm.NWSVel,
                                 vmin=-30,
                                 vmax=30, gatefilter=gatefilter)
                plt.subplot(413)
                display.plot_ppi('sim_velocity', 
                                 sweep=0,
	                         cmap=pyart.graph.cm.NWSVel,
                                 vmin=-30,
                                 vmax=30) 
                plt.savefig(out_path + out_file)
                plt.subplot(414)
                display.plot_ppi('velocity_texture', 
                                 sweep=0,
	                         cmap=pyart.graph.cm.NWSVel,
                                 vmin=0,
                                 vmax=10) 
                plt.savefig(out_path + out_file) 
                plt.close() 
            except:
                import sys
                print('Skipping corrupt time' +
                      year_str + 
                      '-' +
                      month_str + 
                      ' ' + 
                      hour_str + 
                      ':' +
                      minute_str)
                print('Exception: ' + str(sys.exc_info()[0]) + str(sys.exc_info()[1]))
                print(str(sys.exc_info()[2].tb_lineno))


times,dates = time_procedures.get_radar_times_berr_cfradial(start_year, 
                                                            start_month,
                                                            start_day,
                                                            start_hour, 
                                                            start_minute,
                                                            end_year, 
                                                            end_month,
                                                            end_day,
                                                            end_hour, 
                                                            end_minute,
                                                            )
print(dates)

# Go through all of the scans
#for rad_time in times:
#    display_time(rad_time)

if __name__ == '__main__':
    serial = 0

    if(serial == 0):
        # Get iPython cluster
        #Map the code and input to all workers
        #client = Client(scheduler_file=)
        result = db.from_sequence(dates)
        result.map(display_time).compute()

        #Reduce the result to get a list of output
        tt = time.time() - t1
        print(tt)
        print(tt/len(times))
    else:
        for rad_date in dates:
            display_time(rad_date)


    
