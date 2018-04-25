import matplotlib
matplotlib.use('Agg')
import pyart
from netCDF4 import Dataset
import xarray
import numpy as np
from datetime import datetime, timedelta
from copy import deepcopy
import glob
import math
import dask.array as da
import time
import sys
import os
from scipy import interpolate, ndimage
from distributed import Client, LocalCluster, wait


exclude_fields = ['temperature', 'height', 'unfolded_differential_phase', 'specific_attenuation_reflectivity', 'specific_attenuation_differential_reflectivity', 'radar_echo_classification', 'radar_estimated_rain_rate', 'D0', 'NW', 'velocity', 'region_dealias_velocity', 'total_power', 'cross_correlation_ratio', 'differential_reflectivity', 'corrected_differential_reflectivity', 'differential_phase', 'corrected_differential_phase', 'corrected_specific_differential_phase', 'spectrum_width', 'signal_to_noise_ratio', 'ROI']


def get_visst_from_time(cur_time):
    year_str = "%04d" % cur_time.year
    day_str = "%02d" % cur_time.day
    month_str = "%02d" % cur_time.month
    data_list = glob.glob(visst_data_path +
                          'twpvisstpx04*' +
                          year_str +
                          month_str +
                          day_str +
                          '*.cdf')
    if(data_list):
        return Dataset(data_list[0])
    else:
        return []

# Get a Radar object given a time period in the CPOL dataset
def get_grid_from_dda(time):
    year_str = "%04d" % time.year
    month_str = "%02d" % time.month
    day_str = "%02d" % time.day
    hour_str = "%02d" % time.hour
    minute_str = "%02d" % time.minute
    second_str = "%02d" % time.second
    file_name_str = (data_path +
                    'cf_compliant_grid' +
                     year_str +
                     month_str +
                     day_str +
                     hour_str +
                     minute_str + '.nc')
    
    radar = pyart.io.read_grid(file_name_str)
    return radar

def dms_to_decimal(deg, minutes, seconds):
    return deg+minutes/60+seconds/3600

# Convert seconds to midnight to a string format
def seconds_to_midnight_to_string(time_secs_after_midnight):
    hours = math.floor(time_secs_after_midnight/3600)
    minutes = math.floor((time_secs_after_midnight - hours*3600)/60)
    temp = datetime.time(int(hours), int(minutes), )
    return temp.strftime('%H%M%S')

def seconds_to_midnight_to_hm(time_secs_after_midnight):
    hours = math.floor(time_secs_after_midnight/3600)
    minutes = math.floor((time_secs_after_midnight - hours*3600)/60)
    return hours, minutes

# get_grid_times_cpol
#     start_year = Start year of animation
#     start_month = Start month of animation
#     start_day = Start day of animation
#     start_hour = Start hour of animation
#     end_year = End year of animation
#     end_month = End month of animation
#     end_day = End day of animation
#     end_minute = End minute of animation
#     minute_interval = Interval in minutes between scans (default is 5)
# This procedure acquires an array of Grid classes between start_time and end_time  
def get_grid_times_cpol(start_year, start_month, start_day,
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
    date_list_final = []
    for i in days:
        year_str = "%04d" % cur_time.year
        day_str = "%02d" % cur_time.day
        month_str = "%02d" % cur_time.month       
        dir_str = year_str + '/' + year_str + month_str + day_str + '/'
        format_str = (cpol_grid_data_path +
                      dir_str + 
                      'CPOL_' +
                      year_str +
                      month_str + 
                      day_str + 
                      '*' +
                      '.nc')
       
        print('Looking for files with format ' + format_str)
          
        data_list = glob.glob(format_str)
        if(len(data_list) > 0):
            day = datetime(cur_time.year, cur_time.month, cur_time.day, 0, 0, 1)
            date_list_final.append(day)

        for j in range(0, len(data_list)):
            file_list.append(data_list[j])
        cur_time = cur_time + timedelta(days=1)
    
    # Parse all of the dates and time in the interval and add them to the time list
    past_time = []
    for file_name in file_list:
        date_str = file_name[-28:-15]
        year_str = date_str[0:4]
        month_str = date_str[4:6]
        day_str = date_str[6:8]
        hour_str = date_str[9:11]
        minute_str = date_str[11:13]
        cur_time = datetime(int(year_str),
                            int(month_str),
                            int(day_str),
                            int(hour_str),
                            int(minute_str),
                            )
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


# Get a Radar object given a time period in the CPOL dataset
def get_grid_from_cpol(time):
    from datetime import timedelta, datetime
    year_str = "%04d" % time.year
    month_str = "%02d" % time.month
    day_str = "%02d" % time.day
    hour_str = "%02d" % time.hour
    minute_str = "%02d" % time.minute
    
    file_name_str = (cpol_grid_data_path + 
                     year_str +
                     '/' + 
                     year_str +
                     month_str +
                     day_str +
                     '/' +
                     'CPOL_' +
                     year_str +
                     month_str +
                     day_str + '_' +
                     hour_str +
                     minute_str +
                     '_GRIDS_2500m.nc')
    print(file_name_str)
    radar = pyart.io.read_grid(file_name_str, exclude_fields=exclude_fields)
    return radar


def get_echotop_heights(cur_time):
    # First, get VISST Tb 
    echo_top_temps_cpol = []

    try:
        pyart_grid = get_grid_from_cpol(cur_time)
    except:
        print('Py-ART grid not found!')
        return []
    try:
        texture = pyart_grid.fields['velocity_texture']['data']
        z = pyart_grid.fields['reflectivity']['data']
        grid_z = pyart_grid.point_z['data']
        grid_y = pyart_grid.point_y['data']
        grid_x = pyart_grid.point_x['data']
    except:
        return []
    array_shape = texture.shape
    echo_top = np.zeros((array_shape[1],array_shape[2]))
    z_values, y_values, x_values = np.meshgrid(range(0,array_shape[0]),
                                               range(0,array_shape[1]),
                                               range(0,array_shape[2]),
                                               indexing='ij')
    labels = y_values*array_shape[2] + x_values
    in_cloud = np.ma.masked_where(np.logical_or(z.mask == True, texture > 3), texture)
    in_cloud[~in_cloud.mask] = labels[~in_cloud.mask]
    echo_top = ndimage.measurements.maximum(grid_z,
                                            labels=in_cloud,
                                            index=in_cloud)
    echo_top = echo_top[0,:,:]
                
    # Exclude values < 15 km from radar
    dist_from_radar = np.sqrt(np.square(grid_x[0]) + np.square(grid_y[0]))                  
    echo_top = np.ma.masked_where(np.logical_or(echo_top == 0, dist_from_radar < 15000), echo_top)
    Lon_cpol = pyart_grid.point_longitude['data'][0]
    Lat_cpol = pyart_grid.point_latitude['data'][0]

    ds = xarray.Dataset({'ETH': (['y', 'x'], echo_top)},
                             coords={'lon': (['y', 'x'], Lon_cpol),
                                     'lat': (['y', 'x'], Lat_cpol),
                                     'time': cur_time,
                                     'reference_time': cur_time},
                         attrs={'units': 'm', 'long_name': ('CPOL echo top' +
                                                            ' height')})
    the_path = (echo_tops_path + "%04d" % cur_time.year + '/' + "%04d" % cur_time.year +
        "%02d" % cur_time.month + "%02d" % cur_time.day + '/')

    if(not os.path.isdir(the_path)):
        os.makedirs(the_path)
    ds.to_netcdf(path=(the_path + 
                           'echo_tops_' +
		           cur_time.strftime('%Y%m%d%H%M')
                           + '.cdf'), mode='w')
    print(str(cur_time) + ' sucessfully processed!')


if __name__ == '__main__':
    # Input the range of dates and time wanted for the collection of images
    start_year = int(sys.argv[1])
    start_day = int(sys.argv[3])
    start_month = int(sys.argv[2])
    start_hour = 0
    start_minute = 3
    start_second = 0

    end_year = int(sys.argv[4])
    end_month = int(sys.argv[5])
    end_day = int(sys.argv[6])
    end_hour = 3
    end_minute = 1
    end_second = 0

    # Start a cluster with x workers
    cluster = LocalCluster(n_workers=36)
    client = Client(cluster)

    data_path = '/lcrc/group/earthscience/rjackson/multidop_grids/'
    visst_data_path = '/lcrc/group/earthscience/rjackson/visst/'
    echo_tops_path = '/lcrc/group/earthscience/rjackson/echo_tops/echo_tops_150km_2500m/'
    cpol_grid_data_path = '/lcrc/group/earthscience/radar/CPOL_level_1b/GRIDDED/GRID_150km_2500m/'        
    # Get the multidop grid times
    times = get_grid_times_cpol(start_year, start_month, 
	                        start_day, start_hour, 
	                        start_minute, end_year,
	                        end_month, end_day, 
	                        end_hour, end_minute)

    # Calculate PDF
    num_levels = 1
    print('Doing parallel grid loading...')
    import time
    tbs = []
    num_times = len(times)
    print('Processing ' + str(num_times) + ' grids.')
    get_echotop_heights(times[0])
    result = client.map(get_echotop_heights, times)
    wait(result)
    


