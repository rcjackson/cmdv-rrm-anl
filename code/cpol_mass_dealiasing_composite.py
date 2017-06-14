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
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from csu_radartools import csu_kdp
from numba import jit
import sys

# File paths
berr_data_file_path = '/lcrc/group/earthscience/radar/stage/radar_disk_two/berr_rapic/'
data_path_cpol = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_rapic/'
out_file_path = '/lcrc/group/earthscience/rjackson/quicklook_plots/cpol/'

## Berrima - 2009-2011 (new format), 2005-2005 (old format)
start_year = int(sys.argv[1])
start_month = int(sys.argv[2])
start_day = int(sys.argv[3])
start_hour = 1
start_minute = 1

end_year = int(sys.argv[4])
end_month = int(sys.argv[5])
end_day = int(sys.argv[6])
end_hour = 1
end_minute = 2

def display_time(rad_date):
    import pyart
    import matplotlib
    import sys
    sys.path.append('/home/rjackson/cmdv-rrm-anl/code/')
    import time_procedures
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import os
    from datetime import timedelta
    from scipy.stats import gaussian_kde
    from scipy.signal import argrelextrema
    from numba import jit

    # Get a Radar object given a time period in the CPOL dataset
    data_path_cpol = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_rapic/'
    out_file_path = '/lcrc/group/earthscience/rjackson/quicklook_plots/cpol/'
    out_data_path = '/lcrc/group/earthscience/rjackson/cpol/'

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
    def extract_unmasked_data(radar, field, bad=-32768):
        """Simplify getting unmasked radar fields from Py-ART"""
        return deepcopy(radar.fields[field]['data'].filled(fill_value=bad))

    @jit(nopython=True, cache=True)
    def get_fold_position(the_phidp):
 
        tmp = the_phidp
        rth_pos = np.zeros((tmp.shape[0]), dtype=np.int32)  
        for j in range(tmp.shape[0]):
            max_phidp = -32768  
            for i in range(70, tmp.shape[1]):
                if the_phidp[j,i] < max_phidp-30.0 and the_phidp[j,i] > -1000.0:
                    rth_pos[j] = i                    
                    break
                if(the_phidp[j,i] > max_phidp):
                    max_phidp = the_phidp[j,i]
        return rth_pos


    @jit(nopython=True, cache=True)
    def unfold_phidp(the_phidp, rth_position):
        tmp = the_phidp
        for j in range(len(rth_position)):
            i = rth_position[j]
            if i == 0:
                continue
            else:
                tmp[j, i:] += 180
        return tmp

    @jit(cache=True)
    def refold_vdop(vdop_art, v_nyq_vel, rth_position):
        tmp = vdop_art
        for j in range(len(rth_position)):
            i = rth_position[j]
            if i == 0:
                continue
            else:
                tmp[j, i:] += v_nyq_vel
    
        pos = (vdop_art > v_nyq_vel)
        tmp[pos] = tmp[pos] - 2*v_nyq_vel
        return tmp
        
    one_day_later = rad_date+timedelta(days=1)
    times, dates = time_procedures.get_radar_times_cpol_cfradial(rad_date.year, 
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
    rerun_multidop_list = open('/home/rjackson/grids_to_regenerate/folds' + 
                               str(rad_date.year) + 
                               str(rad_date.month) + 
                               str(rad_date.day), mode='w+')
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

        try:  
            radar = time_procedures.get_radar_from_cpol_cfradial(rad_time)
            if(rad_time.year > 2007):
                ref_field = 'Refl'
                vel_field = 'Vel'
                ref_threshold = 0
                rhohv_field = 'RHOHV'
                phidp_field = 'PHIDP'
            else:
                ref_field = 'reflectivity'
                vel_field = 'velocity'
                ref_threshold = 5
                rhohv_field = 'cross_correlation_ratio'
                phidp_field = 'differential_phase'
            if(radar.nsweeps == 1):
                raise Exception('Radar only has one sweep!')
            if(not 'last_Radar' in locals()):
                    last_Radar = radar

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
            steps = np.floor(len(u)/50)
            wind_profile = pyart.core.HorizontalWindProfile.from_u_and_v(alt[0::steps],
                                                                         u[0::steps],
                                                                         v[0::steps])
      
            ## 4DD expects speed, direction but HorizontalWindProfile outputs u_wind, v_wind
            wind_profile.u = wind_profile.u_wind
            wind_profile.v = wind_profile.v_wind
            nyq = radar.instrument_parameters['nyquist_velocity']['data'][0]   
            print(wind_profile)
            # Dealias velocities
            gatefilter = pyart.correct.GateFilter(radar)         
            # Region based hangs before 2006 hangs on noise w/Z < 10 dBZ
            gatefilter.exclude_below(ref_field, ref_threshold)
            try:
                gatefilter.exclude_below(rhohv_field, 0.6)
            except KeyError: 
                print('No rho HV data! Disabling Rho HV mask!')
            gatefilter.exclude_masked(vel_field)
            gatefilter.exclude_invalid(vel_field)
            gatefilter.exclude_masked(ref_field)
            gatefilter.exclude_invalid(ref_field)
            gatefilter.exclude_above(ref_field, 80)
            gatefilter.exclude_below(vel_field, -75)
            gatefilter.exclude_above(vel_field, 75)
            # Look for phiDP folds and correct velocities
            try:
                dz = radar.fields[ref_field]['data']
                dp = radar.fields[phidp_field]['data']
                r, azi = radar.range['data'], radar.azimuth['data']
                rng2d, az2d = np.meshgrid(radar.range['data'], radar.azimuth['data'])
                kdN, fdN, sdN = csu_kdp.calc_kdp_bringi(dp=dp, dz=dz, rng=rng2d/1000.0, thsd=24, gs=250.0, window=4)

                fdN = np.ma.MaskedArray(fdN)
                fdN[fdN == -32768] = np.ma.masked
                  
                phidp_fold = deepcopy(fdN)
                phidp_fold[phidp_fold.mask == True] = np.nan
                rth = get_fold_position(phidp_fold)
               
                vdop_art = deepcopy(radar.fields[vel_field]['data'])
                v_nyq_vel = radar.instrument_parameters['nyquist_velocity']['data'][0]
                
                vdop_refolded = refold_vdop(vdop_art, v_nyq_vel, rth)
                radar.add_field_like(vel_field, 'vdop_phidp_fold_corrected', 
                                     vdop_refolded, replace_existing=True)
                radar.add_field_like(phidp_field, 'phidp_bringi', 
                                     phidp_fold, replace_existing=True)  
                vel_field = 'vdop_phidp_fold_corrected'
                 
                # Tell me if fold is within 60 km of CPOL and write to file (for multidop reprocessing purposes)
                is_fold = 0
                for fold_starts in radar.range['data'][rth]:
                    if(fold_starts > 13500 and fold_starts < 80000 and is_fold == 0):
                        print('Fold in DD domain at ' + str(fold_starts/1e3) + ' km')
                        print(rad_time.strftime('%m/%j/%Y %H:%M:%S'))
                        rerun_multidop_list.write(hour_str + ':' + minute_str + '/n')
                        is_fold = 1
                gatefilter.exclude_masked('phidp_bringi')
            except:
                import sys
                exc_type, exc_obj, exc_tb = sys.exc_info()
                print('Cant run phiDP unfolding!')
                print('Exception: ' + str(sys.exc_info()[0]) 
                                    + str(sys.exc_info()[1])
                                    + str(exc_tb.tb_lineno))
 

            print('Running 4DD...')
            vels = pyart.correct.dealias._create_rsl_volume(radar, 
                                                            vel_field, 
                                                            0, 
                                                            -9999.0, 
                                                            excluded=None)

            # Discontinuties in azimuthal angles due to corrupt data make 4DD segfault
            for i in range(0,radar.nsweeps-1):
                sweep = vels.get_sweep(i)
                ray0 = sweep.get_ray(0)
                ray50 = sweep.get_ray(50)
                diff = ray0.azimuth-ray50.azimuth 
                if(diff > 180.0):
                    diff = 360.0 - diff    
                if(abs(diff)/50.0 < 0.5):
                    print('Corrupt azimuthal angle data....skipping file!')
                    raise Exception('Corrupt azimuthal angles!')

            if(last_Radar.nsweeps == radar.nsweeps and not last_Radar == radar):
                try:
                    corrected_velocity_4dd = pyart.correct.dealias_fourdd(radar,
                                                                          vel_field=vel_field,
                                                                          keep_original=False,
                                                                          filt=1,
                                                                          sonde_profile=wind_profile,
                                                                          gatefilter=gatefilter,
                                                                          ) 
                except:
                    corrected_velocity_4dd = pyart.correct.dealias_fourdd(radar,
                                                                          vel_field=vel_field,
                                                                          keep_original=False,
                                                                          filt=1,
                                                                          sonde_profile=wind_profile,
                                                                          gatefilter=gatefilter,
                                                                          ) 
            else:
                corrected_velocity_4dd = pyart.correct.dealias_fourdd(radar,
                                                                      vel_field=vel_field,
                                                                      keep_original=False,
                                                                      filt=1,
                                                                      sonde_profile=wind_profile,
                                                                      gatefilter=gatefilter, 
                                                                      )
            print('Running region based dealiasing...')
            # Calculate result from region based dealiasing
            corrected_velocity_region = pyart.correct.dealias_region_based(radar,
                                                                           vel_field=vel_field,
                                                                           keep_original=False,
                                                                           centered=True,
                                                                           gatefilter=gatefilter,
                                                                           interval_splits=6,
                                                                           skip_between_rays=2000,
                                                                           skip_along_ray=2000,
                                                                           rays_wrap_around=True,
                                                                           valid_min=-75,
                                                                           valid_max=75,
                                                                           nyquist_velocity=nyq)                                              
            corr_vel = corrected_velocity_4dd['data']
            difference = corr_vel/corrected_velocity_region['data']
            corr_vel = np.ma.masked_where(np.logical_or(difference > 1.05, 
                                                        difference < 0.95),
                                          corr_vel)
            corrected_velocity_4dd['data'] = corr_vel
            radar.add_field_like(vel_field, 
                                 'corrected_velocity', 
                                 corrected_velocity_4dd['data'],
	                         replace_existing=True)
 
            # Save to Cf/Radial file
            time_procedures.write_radar_to_cpol_cfradial(radar, rad_time)
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
            plt.figure(figsize=(7,24))
            plt.subplot(311)
            display = pyart.graph.RadarMapDisplay(radar)
            display.plot_ppi(ref_field,
                             gatefilter=gatefilter, 
                             sweep=0, 
                             cmap=pyart.graph.cm.NWSRef,
                             vmin=0, 
                             vmax=70)
            plt.subplot(312)
            display = pyart.graph.RadarMapDisplay(radar)
            display.plot_ppi('corrected_velocity', 
                             gatefilter=gatefilter, 
                             sweep=0,
                             cmap=pyart.graph.cm.NWSVel,
                             vmin=-30,
                             vmax=30)
            plt.savefig(out_path + out_file)
            plt.subplot(313)
            display = pyart.graph.RadarMapDisplay(radar)
            display.plot_ppi('phidp_bringi', 
                             gatefilter=gatefilter,
                             sweep=0,
                             cmap='jet',
                             vmin=-360,
                             vmax=360)
            plt.savefig(out_path + out_file)
            plt.close() 
        except:
            import sys
            exc_type, exc_obj, exc_tb = sys.exc_info()
            print('Skipping corrupt time' +
                  year_str + 
                  '-' +
                  month_str + 
                  ' ' + 
                  hour_str + 
                  ':' +
                  minute_str)
            print('Exception: ' + 
                  str(sys.exc_info()[0]) + 
                  str(sys.exc_info()[1]) +
                  str(exc_tb.tb_lineno))
    rerun_multidop_list.close()    
 
times,dates = time_procedures.get_radar_times_cpol_cfradial(start_year, 
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

serial = 1
if(serial == 0):
    # Go through all of the scans
    #for rad_time in times:
    #    display_time(rad_time)

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

    #Map the code and input to all workers
    result = My_View.map_async(display_time, dates)
    #Reduce the result to get a list of output
    qvps = result.get()
    tt = time.time() - t1
    print(tt)
    print(tt/len(times))
else:
    for date in dates:
        display_time(date)


    
