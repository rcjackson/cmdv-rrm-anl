import matplotlib
matplotlib.use('Agg')
import pyart
from matplotlib import pyplot as plt
import numpy as np
import glob
import os
from copy import deepcopy
from ipyparallel import Client
from time import sleep
import time
import time_procedures
import sys

# File paths
berr_data_file_path = '/lcrc/group/earthscience/radar/stage/radar_disk_two/berr_rapic/'
data_path_cpol = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_rapic/'
out_file_path = '/lcrc/group/earthscience/rjackson/quicklook_plots/berr/'

## Berrima - 2009-2011 (new format), 2005-2005 (old format)
start_year = 2005
start_month = 12
start_day = 24
start_hour = 0
start_minute = 50 

end_year = 2005
end_month = 12
end_day = 25
end_hour = 0
end_minute = 2

def display_time(rad_time):
    import pyart
    import matplotlib
    import sys
    sys.path.append('/home/rjackson/cmdv-rrm-anl/code/')
    import time_procedures
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import os
    from datetime import timedelta

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
    
    
 
    year_str = "%04d" % rad_time.year
    month_str = "%02d" % rad_time.month
    day_str = "%02d" % rad_time.day
    hour_str = "%02d" % rad_time.hour
    minute_str = "%02d" % rad_time.minute
    second_str = "%02d" % rad_time.second
  
    try:
        radar = time_procedures.get_radar_from_berr_cfradial(rad_time)
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
        steps = np.floor(len(u)/30)
        wind_profile = pyart.core.HorizontalWindProfile.from_u_and_v(alt[0::steps],
                                                                     u[0::steps],
                                                                     v[0::steps])
              
        ## 4DD expects speed, direction but HorizontalWindProfile outputs u_wind, v_wind
        wind_profile.u = wind_profile.u_wind
        wind_profile.v = wind_profile.v_wind
         
        # Dealias velocities
        gatefilter = pyart.correct.GateFilter(radar)
        gatefilter.exclude_below(ref_field, 5)
        gatefilter.exclude_invalid(vel_field)
        gatefilter.exclude_masked(vel_field)
        gatefilter.exclude_invalid(ref_field) 
        gatefilter = pyart.correct.despeckle_field(radar,
                                                   vel_field, 
                                                   gatefilter=gatefilter,
                                                   size=50)     
        radar.add_field('sim_velocity',
                        pyart.util.simulated_vel_from_profile(radar, 
                                                              wind_profile),
                        replace_existing = True)
        nyq = radar.instrument_parameters['nyquist_velocity']['data'][1]
        
        corrected_velocity_4dd = pyart.correct.dealias_region_based(radar,
                                                                    vel_field=vel_field,
                                                                    keep_original=False,
                                                                    centered=True,
                                                                    gatefilter=gatefilter,
                                                                    interval_splits=6,
                                                                    skip_between_rays=2200,
                                                                    skip_along_ray=2200,
                                                                    rays_wrap_around=True,
                                                                    valid_min=-75,
                                                                    valid_max=75,
                                                                    nyquist_velocity=nyq)
              
        # Filter out regions based on deviation from sounding field to remove missed folds                                                         
        corr_vel = corrected_velocity_4dd['data']
        sim_velocity = radar.fields['sim_velocity']['data']
        diff = corr_vel - sim_velocity
        diff = diff/(nyq)                   
        radar.add_field_like(vel_field, 
                             'corrected_velocity', 
  	                     corrected_velocity_4dd['data'],
	                     replace_existing=True)
 
        # Calculate gradient of field
        gradient = pyart.config.get_metadata('velocity')
        gradients = np.ma.array(np.gradient(radar.fields['corrected_velocity']['data']))
        gradients = np.ma.masked_where(gradients < -31000,gradients)
        gradients = gradients/nyq
        gradient['data'] = gradients[0]
        gradient['standard_name'] = 'gradient_of_corrected_velocity_wrt_azimuth'
        gradient['units'] = 'meters per second per gate (divided by Vn)'
        radar.add_field('gradient_wrt_angle',
                        gradient,
                        replace_existing=True)

        gradient = pyart.config.get_metadata('velocity')
        gradient['data'] = gradients[1]
        gradient['standard_name'] = 'gradient_of_corrected_velocity_wrt_range'
        gradient['units'] = 'meters per second per gate (divided by Vn)'
        radar.add_field('gradient_wrt_range',
                        gradient,
                        replace_existing=True)

        # Calculate difference from simulated velocity  
        radar.add_field_like('sim_velocity', 
                             'velocity_diff', 
                             diff, 
                             replace_existing=True) 
        
        # Adjust sweeps to match reference velocity
        ref_vdata = sim_velocity
        corr_vel = corrected_velocity_4dd['data']
        nyquist_interval = float(2*nyq)
        for nsweep, sweep_slice in enumerate(radar.iter_slice()):                                                                   
            sref = ref_vdata[sweep_slice]
            scorr = corr_vel[sweep_slice]
            mean_diff = (sref - scorr).mean()
            if(mean_diff > -100):
                global_fold = round(mean_diff / nyquist_interval)
                if global_fold != 0:
                    corr_vel[sweep_slice] += global_fold * nyquist_interval
            
        # Filter by devation from sounding  
        fifth_percentile = np.percentile(diff[diff > -100], 10) 
        ninety_fifth_percentile = np.percentile(diff[diff > -100], 90) 
        if(ninety_fifth_percentile > 2.5):
            ninety_fifth_percentile = 2.5
        if(fifth_percentile < -1.0):
            fifth_percentile = -1.0
        inds2 = np.logical_or(diff < fifth_percentile, diff > ninety_fifth_percentile)
        inds1 = np.logical_and(corr_vel < -nyq/2, sim_velocity > nyq/2)
        inds2 = np.logical_and(corr_vel > nyq/2, sim_velocity < -nyq/2)
        corr_vel = np.ma.masked_where(np.logical_or(inds1,inds2), corr_vel)        
        corrected_velocity_4dd['data'] = corr_vel
        radar.add_field_like(vel_field, 
                             'corrected_velocity', 
  	                     corrected_velocity_4dd['data'],
	                     replace_existing=True)
	gatefilter = pyart.correct.despeckle_field(radar,
                                                   'corrected_velocity', 
                                                   size=50)         
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
            try:
                os.makedirs(out_path)
            except:
                print('Not making directory')
        out_file = hour_str + minute_str + '.png'
        plt.figure(figsize=(7,14))
        plt.subplot(211)
        display = pyart.graph.RadarMapDisplay(radar)
        display.plot_ppi(ref_field, 
                         sweep=0, 
                         cmap=pyart.graph.cm.NWSRef,
                         vmin=0, 
                         vmax=70)
        plt.subplot(212)
        display = pyart.graph.RadarMapDisplay(radar)
        display.plot_ppi('corrected_velocity', 
                         sweep=0,
	                 cmap=pyart.graph.cm.NWSVel,
                         vmin=-30, gatefilter=gatefilter,
                         vmax=30)
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

serial = 1

if(serial == 0):
    # Get iPython cluster
    state = 0
    while state == 0:
        try:
            My_Cluster = Client()
            My_View = My_Cluster[:]
            state = 1
        except:
            state = 0
            print('Cluster not ready for me')
            sleep(10)

    #Turn off blocking so all engines can work async
    My_View.block = False

    #on all engines do an import of Py-ART
    My_View.execute('import matplotlib')
    My_View.execute('matplotlib.use("agg")')
    My_View.execute('import pyart')
    My_View.execute('import numpy as np')
    My_View.execute('import time_procedures')
    t1 = time.time()

    #for rad_date in dates:
    #    display_time(rad_date)

    #Map the code and input to all workers
    result = My_View.map_async(display_time, times)

    #Reduce the result to get a list of output
    qvps = result.get()
    tt = time.time() - t1
    print(tt)
    print(tt/len(times))
else:
    for rad_date in times:
        display_time(rad_date)


    
