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

# File paths
berr_data_file_path = '/lcrc/group/earthscience/radar/stage/radar_disk_two/berr_rapic/'
data_path_cpol = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_rapic/'
out_file_path = '/lcrc/group/earthscience/rjackson/quicklook_plots/cpol/'

## Berrima - 2009-2011 (new format), 2005-2005 (old format)
start_year = 2005
start_month = 12
start_day = 24
start_hour = 0
start_minute = 1

end_year = 2005
end_month = 12
end_day = 25
end_hour = 0
end_minute = 2

# get_radar_times_cpol
#     start_year = Start year of animation
#     start_month = Start month of animation
#     start_day = Start day of animation
#     start_hour = Start hour of animation
#     end_year = End year of animation
#     end_month = End month of animation
#     end_day = End day of animation
#     end_minute = End minute of animation
#     minute_interval = Interval in minutes between scans (default is 5)
# This procedure acquires an array of Radar classes between start_time and end_time  
def get_radar_times_cpol(start_year, start_month, start_day,
                         start_hour, start_minute, end_year,
                         end_month, end_day, end_hour, 
                         end_minute, minute_interval=5):

    from datetime import timedelta, datetime
    start_time = datetime(start_year,
                          start_month,
                          start_day,
                          start_hour,
                          start_minute,
                          )
    end_time = datetime(end_year,
                        end_month,
                        end_day,
                        end_hour,
                        end_minute,
                        )  

    deltatime = end_time - start_time
    

    if(deltatime.seconds > 0 or deltatime.minute > 0):
        no_days = deltatime.days + 1
    else:
        no_days = deltatime.days
    
    if(start_day != end_day):
        no_days = no_days + 1
        
    days = range(0, no_days)
    print('We are about to load grid files for ' + str(no_days) + ' days')
    
    

    # Find the list of files for each day
    cur_time = start_time
 
    file_list = []
    time_list = []
    for i in days:
        year_str = "%04d" % cur_time.year
        day_str = "%02d" % cur_time.day
        month_str = "%02d" % cur_time.month
        if(cur_time.year == 2009 or (cur_time.year == 2010 and
                                     cur_time.month < 6 )):
            dir_str = 'cpol_0910/rapic/'
        else:
            dir_str = 'cpol_1011/rapic/'

        format_str = (data_path_cpol +
                      dir_str +
                      year_str +
                      month_str + 
                      day_str + 
                      '*Gunn_Pt' +
                      '.rapic')
    
    
        print('Looking for files with format ' + format_str)
          
        data_list = glob.glob(format_str)
        
        for j in range(0, len(data_list)):
            file_list.append(data_list[j])
        cur_time = cur_time + timedelta(days=1)
    
    # Parse all of the dates and time in the interval and add them to the time list
    past_time = []
    for file_name in file_list:
        date_str = file_name[-25:-11]
        year_str = date_str[0:4]
        month_str = date_str[4:6]
        day_str = date_str[6:8]
        hour_str = date_str[8:10]
        minute_str = date_str[10:12]
        
        
        
        
        cur_time = datetime(int(year_str),
                            int(month_str),
                            int(day_str),
                            int(hour_str),
                            int(minute_str),
                            0)
        time_list.append(cur_time)
    
    # Sort time list and make sure time are at least xx min apart
    time_list.sort()
    time_list_sorted = deepcopy(time_list)
   
    time_list_final = []
    past_time = []
    
    
    for times in time_list_sorted: 
        
        cur_time = times  
        
        if(past_time == []):
            past_time = cur_time
            
        if(cur_time - past_time >= timedelta(minutes=minute_interval)
           and cur_time >= start_time and cur_time <= end_time):
            
            time_list_final.append(cur_time)
            past_time = cur_time
                   
    return time_list_final



def display_time(rad_time):
    import pyart
    import matplotlib
    import os
    os.chdir('/home/rjackson/cmdv-rrm-anl/code/')
    import time_procedures
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import os
    from datetime import timedelta
    from scipy import ndimage
    

    # Get a Radar object given a time period in the CPOL dataset
    data_path_cpol = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_rapic/'
    out_file_path = '/lcrc/group/earthscience/rjackson/quicklook_plots/cpol/'
    out_data_path = '/lcrc/group/earthscience/rjackson/cpol/'

    # CPOL in lassen or rapic?
    cpol_format = 0    # 0 = lassen, 1 = rapic

    year_str = "%04d" % rad_time.year
    month_str = "%02d" % rad_time.month
    day_str = "%02d" % rad_time.day
    hour_str = "%02d" % rad_time.hour
    minute_str = "%02d" % rad_time.minute
    second_str = "%02d" % rad_time.second
    
   
    try:
        radar = time_procedures.get_radar_from_cpol_cfradial(rad_time)
        if(rad_time.year > 2007):
            ref_field = 'Refl'
            vel_field = 'Vel'
            rhohv_field = 'RHOHV'
        else:
            ref_field = 'reflectivity'
            vel_field = 'velocity'
            rhohv_field = 'cross_correlation_ratio'

        if(radar.nsweeps == 1):
            return
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
        # Filter clutter and noise from velocities

        nyq_Gunn = radar.instrument_parameters['nyquist_velocity']['data'][0]
	    
        gatefilter = pyart.correct.GateFilter(radar)
        gatefilter.exclude_below(ref_field, 10)
        gatefilter.exclude_below(rhohv_field, 0.5)
        gatefilter.exclude_invalid(vel_field)
        gatefilter.exclude_masked(vel_field)
        gatefilter.exclude_invalid(ref_field)
        gatefilter = pyart.correct.despeckle_field(vel_field, 
                                                   gatefilter=gatefilter,
                                                   size=10) 
        radar.add_field('sim_velocity',
                        pyart.util.simulated_vel_from_profile(radar, 
                                                              wind_profile),
                        replace_existing = True)
	    
        corrected_velocity_4dd = pyart.correct.dealias_region_based(radar,
                                                                    vel_field=vel_field,
                                                                    keep_original=False,
	                                                            centered=True,
	                                                            interval_splits=6,
                                                                    gatefilter=gatefilter,
	                                                            skip_between_rays=2000,
	                                                            skip_along_ray=2000,
	                                                            rays_wrap_around=True,
	                                                            valid_min=-75,
	                                                            valid_max=75)
        
        print('Dealiasing done!')
        # Filter out regions based on deviation from sounding field to remove missed folds                                                         
        corr_vel = corrected_velocity_4dd['data']
        sim_velocity = radar.fields['sim_velocity']['data']
        diff = corr_vel - radar.fields['sim_velocity']['data']
        diff = diff/(radar.instrument_parameters['nyquist_velocity']['data'][1])                   
        radar.add_field_like(vel_field, 
                             'corrected_velocity', 
  	                     corrected_velocity_4dd['data'],
	                     replace_existing=True)
 
        # Calculate gradient of field
        gradient = pyart.config.get_metadata('velocity')
        gradients = np.ma.array(np.gradient(radar.fields['corrected_velocity']['data']))
        gradients = np.ma.masked_where(gradients < -31000,gradients)
        gradients = gradients/(radar.instrument_parameters['nyquist_velocity']['data'][1])
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

        # Adjust sweeps to match reference velocity
        ref_vdata = sim_velocity
        corr_vel = corrected_velocity_4dd['data']
        nyquist_interval = float(2*radar.instrument_parameters['nyquist_velocity']['data'][1])
        for nsweep, sweep_slice in enumerate(radar.iter_slice()):                                                                   
            sref = ref_vdata[sweep_slice]
            scorr = corr_vel[sweep_slice]
            mean_diff = (sref - scorr).mean()
            if(mean_diff > -100):
                global_fold = round(mean_diff / nyquist_interval)
                if global_fold != 0:
                    corr_vel[sweep_slice] += global_fold * nyquist_interval
       
        # Calculate difference from simulated velocity
        diff = radar.fields['corrected_velocity']['data'] - radar.fields['sim_velocity']['data']
        diff = diff/(radar.instrument_parameters['nyquist_velocity']['data'][1])     
        radar.add_field_like('sim_velocity', 
                             'velocity_diff', 
                             diff, 
                             replace_existing=True)    
      
        # Filter by gradient   
        corr_vel = np.ma.masked_where(np.logical_or(np.logical_or(gradients[0] > 0.3, 
                                                                  gradients[0] < -0.3),
                                                    np.logical_or(diff > 2.0, 
                                                                  diff < -0.9)),
                                      corr_vel)
        corrected_velocity_4dd['data'] = corr_vel
        radar.add_field_like(vel_field, 
                             'corrected_velocity', 
  	                     corrected_velocity_4dd['data'],
	                     replace_existing=True)
	print('Filter!')    
        # Save to Cf/Radial file
        time_procedures.write_radar_to_cpol_cfradial(radar, rad_time)
        
        out_path = (out_file_path +
                    '/' +
                    year_str +
                    '/' +
                    month_str +
                    '/' +
                    day_str +
                    '/')
        if(not os.path.exists(out_path)):
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
                         vmin=-30,
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
My_View.execute('from time_procedures import get_radar_from_cpol_rapic')
t1 = time.time()

#for timer in times:
#    display_time(timer)

#Map the code and input to all workers
result = My_View.map_async(display_time, times)

#Reduce the result to get a list of output
qvps = result.get()
tt = time.time() - t1
print(tt)
print(tt/len(times))


    
