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
from scipy import interpolate, ndimage
import time_procedures
from distributed import Client, LocalCluster
from dask import delayed, compute

# Input the range of dates and time wanted for the collection of images
start_year = int(sys.argv[2])
start_day = int(sys.argv[4])
start_month = int(sys.argv[3])
start_hour = 0
start_minute = 1
start_second = 0

end_year = int(sys.argv[5])
end_month = int(sys.argv[6])
end_day = int(sys.argv[7])
end_hour = 3
end_minute = 1
end_second = 0

# Start a cluster with x workers
cluster = LocalCluster(n_workers=int(sys.argv[1]))
client = Client(cluster)

data_path = '/lcrc/group/earthscience/rjackson/multidop_grids/ddop/'
visst_data_path = '/lcrc/group/earthscience/rjackson/visst/'
echo_tops_path = '/lcrc/group/earthscience/rjackson/echo_tops/'

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

def get_echotop_heights(cur_time):
    # First, get VISST Tb 
    echo_top_temps_cpol = []

    try:
        pyart_grid = time_procedures.get_grid_from_cpol(cur_time)
    except:
        print('Py-ART grid not found!')
        return []
    texture = pyart_grid.fields['velocity_texture']['data']
    grid_z = pyart_grid.point_z['data']
    grid_y = pyart_grid.point_y['data']
    grid_x = pyart_grid.point_x['data']
    array_shape = texture.shape
    echo_top = np.zeros((array_shape[1],array_shape[2]))
    echo_temp = np.zeros((array_shape[1],array_shape[2]))
    z_values, y_values, x_values = np.meshgrid(range(0,array_shape[0]),
                                               range(0,array_shape[1]),
                                               range(0,array_shape[2]),
                                               indexing='ij')
    labels = y_values*array_shape[2] + x_values
    in_cloud = np.ma.masked_where(texture > 3, texture)
    in_cloud[~in_cloud.mask] = labels[~in_cloud.mask]
    echo_top = ndimage.measurements.maximum(grid_z,
                                            labels=in_cloud,
                                            index=in_cloud)
    echo_top = echo_top[0,:,:]
                    
    # Exclude values < 15 km from radar
    dist_from_radar = np.sqrt(np.square(grid_x[0]) + np.square(grid_y[0]))                  
    echo_top = np.ma.masked_where(dist_from_radar < 15000, echo_top)
    return echo_top
        
# Get the multidop grid times
times = time_procedures.get_grid_times_cpol(start_year, start_month, 
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
hours = []
minutes = []
seconds = []
years = []
days = []
months = []
num_frames = 2000
first_grid = time_procedures.get_grid_from_cpol(times[0])
Lon_cpol = first_grid.point_longitude['data'][0]
Lat_cpol = first_grid.point_latitude['data'][0]
first_array = get_echotop_heights(times[0])
get_heights = delayed(get_echotop_heights)

for cur_time in times:
    tbs_shape = first_array.shape
    years.append(cur_time.year*np.ones(tbs_shape[1]))
    days.append(cur_time.day*np.ones(tbs_shape[1]))
    months.append(cur_time.month*np.ones(tbs_shape[1]))
    hours.append(cur_time.hour*np.ones(tbs_shape[1]))
    minutes.append(cur_time.minute*np.ones(tbs_shape[1]))
    seconds.append(cur_time.second*np.ones(tbs_shape[1]))

years = np.concatenate([arrays for arrays in years])
days = np.concatenate([arrays for arrays in days])
months = np.concatenate([arrays for arrays in months])
hours = np.concatenate([arrays for arrays in hours])
minutes = np.concatenate([arrays for arrays in minutes])
seconds = np.concatenate([arrays for arrays in seconds])

for i in range(0, len(years), num_frames):
    t1 = time.time()
    tbs_temp = [da.from_delayed(get_heights(cur_time), 
                                shape=first_array.shape, 
                                dtype=float) for cur_time in times[i:i+num_frames]]
    tbs_temp = da.stack(tbs_temp, axis=0)
    print(tbs_temp.shape)

    tbs_temp = compute(*tbs_temp)
    tbs_temp = np.stack(tbs_temp)

    ds = xarray.Dataset({'cpol_T': (['time', 'y', 'x'], tbs_temp)},
                         coords={'lon': (['y', 'x'], Lon_cpol),
                                'lat': (['y', 'x'], Lat_cpol),
                                'time': times[i:i+num_frames],
                                'reference_time': times[i]},
                         attrs={'units': 'K', 'long_name': ('CPOL echo top' +  
                                                            ' temperature')})
    print(ds)
    ds.to_netcdf(path=(echo_tops_path +
                       'echo_tops' + 
                       times[i].strftime('%Y%m%j%H%M') 
                       + '.cdf'), mode='w')
    t2 = time.time() - t1
    print('Total time in s: ' + str(t2))
    print('Time per scan = ' + str(t2/num_frames))
        
print(years.shape)
   
print('Writing netCDF file...')


