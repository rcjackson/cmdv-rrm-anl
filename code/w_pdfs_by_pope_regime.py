import matplotlib
matplotlib.use('Agg')
import pyart
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
from copy import deepcopy
import glob
import math
import dask.array as da
from distributed import Client, LocalCluster
from dask import delayed
import time
import sys

# Start a cluster with x workers
cluster = LocalCluster(n_workers=int(sys.argv[1]))
client = Client(cluster)

# Input the range of dates and time wanted for the collection of images
start_year = 2005
start_day = 1
start_month = 11
start_hour = 1
start_minute = 0
start_second = 0

end_year = 2011
end_month = 5
end_day = 1
end_hour = 0
end_minute = 00
end_second = 0

data_path = '/lcrc/group/earthscience/rjackson/multidop_grids/ddop/'


# get_radar_times
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
def get_dda_times(start_year, start_month, start_day,
                  start_hour, start_minute, end_year,
                  end_month, end_day, end_hour, 
                  end_minute, minute_interval=5):

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
        
    days = np.arange(0, no_days, 1)
    print('We are about to load grid files for ' + str(no_days) + ' days')
    
    # Find the list of files for each day
    cur_time = start_time
 
    file_list = []
    time_list = []
    for i in days:
        year_str = "%04d" % cur_time.year
        day_str = "%02d" % cur_time.day
        month_str = "%02d" % cur_time.month
        format_str = (data_path +
                      'cf_compliant_grid' +
                      year_str +
                      month_str +
                      day_str +
                      '*.nc')
    
        print('Looking for files with format ' + format_str)
          
        data_list = glob.glob(format_str)
        
        for j in range(0, len(data_list)):
            file_list.append(data_list[j])
        cur_time = cur_time + timedelta(days=1)
    
    # Parse all of the dates and time in the interval and add them to the time list
    past_time = []
    for file_name in file_list:
        date_str = file_name[-15:-3]
        year_str = date_str[0:4]
        month_str = date_str[4:6]
        day_str = date_str[6:8]
        hour_str = date_str[8:10]
        minute_str = date_str[10:12]
        second_str = '00'
             
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

# Get beam crossing angle between radars
def get_bca(grid):
    berr_origin = [-12960.1,-23091.1]
    x,y = np.meshgrid(grid.x['data'], grid.y['data'])
    a = np.sqrt(np.multiply(x,x)+np.multiply(y,y))
    b = np.sqrt(pow(x-berr_origin[0],2)+pow(y-berr_origin[1],2))
    c = np.sqrt(berr_origin[0]*berr_origin[0]+berr_origin[1]*berr_origin[1])
    theta_1 = np.arccos(x/a)
    theta_2 = np.arccos((x-berr_origin[1])/b)
    return np.arccos((a*a+b*b-c*c)/(2*a*b))

def get_updrafts(time):
    pyart_grid = get_grid_from_dda(time)
    bca = get_bca(pyart_grid)
    w = pyart_grid.fields['upward_air_velocity']['data']
    updraft_depth = np.zeros(w[0].shape)
    z = pyart_grid.fields['reflectivity']['data']
    bca = np.ma.masked_invalid(bca)
    for levels in range(0,num_levels-1):
        outside_dd_lobes = np.logical_or(bca < math.pi/6, 
                                         bca > 5*math.pi/6)
        w[levels] = np.ma.masked_where(np.logical_or(outside_dd_lobes, 
                                                     z[levels] < 1), 
                                       w[levels])
        is_in_updraft = w[levels] > 1
        is_in_updraft_next = w[levels+1] > 1
        both_in_updraft = np.logical_or(np.logical_and(is_in_updraft,
                                                       is_in_updraft_next),
                                        updraft_depth > 10)
        
        add_one = np.where(both_in_updraft)
        set_to_zero = np.where(~both_in_updraft)
        if(len(add_one[0]) > 0):
            updraft_depth[add_one[0], 
                          add_one[1]] = updraft_depth[add_one[0], add_one[1]] + 1 
            updraft_depth[set_to_zero[0], 
                          set_to_zero[1]] = 0  
        
    for levels in range(0,num_levels-1):
        invalid_w = np.logical_or(w[levels] < -99, w[levels] > 99)
        outside_updraft = np.logical_or(updraft_depth < 10, z[levels] < 1)
        outside_updraft_and_lobes = np.logical_or(outside_updraft, 
                                                  outside_dd_lobes)
        w[levels] = np.ma.masked_where(np.logical_or(invalid_w, 
                                                     outside_updraft_and_lobes)
                                       , w[levels])
    w[w.mask == True] = np.nan
        
    # Make new array that is 1 by num_levels by 81 by 111
    ws_temp = np.ma.zeros((1,num_levels,81,111))
    rfs_temp = np.ma.zeros((1,num_levels,81,111))
    ws_temp[0] = w
    rfs_temp[0] = z 
    return w

# Plot the radars from given time.

times = get_dda_times(start_year, start_month, start_day,
                      start_hour, start_minute, end_year,
                      end_month, end_day, end_hour, 
                      end_minute, minute_interval=0)

in_netcdf = Dataset('/lcrc/group/earthscience/rjackson/data/Pope_regime.cdf', 
                    mode='r')            
year = in_netcdf.variables['year'][:]
month = in_netcdf.variables['month'][:]
day = in_netcdf.variables['day'][:]
groups = in_netcdf.variables['groups'][:]

popedates = []
for i in range(0,len(day)):
    popedates.append(datetime(year=int(year[i]),
                              month=int(month[i]),
                              day=int(day[i])))

# Since grids are uniform, calculate beam crossing angle for first grid and
# apply to all
first_grid = get_grid_from_dda(times[0])
bca = get_bca(first_grid) 
num_levels = 40
z_levels = np.arange(0.5,0.5*(num_levels+1),0.5)*1000
count = 0
pope_regime = int(sys.argv[2])

# Filter out data not in Pope regime
pope_times = []
for time in times:
    # Look for date in Pope regime data
    cur_date = datetime(year=time.year, month=time.month, day=time.day)
    inds = np.where([day <= cur_date for day in popedates])
    pope_index = inds[0][-1]
    print((popedates[pope_index], time))
    if(groups[pope_index] == pope_regime):
        pope_times.append(time)

in_netcdf.close()

# Get delayed structure to load files in parallel
get_file = delayed(get_updrafts)
ws = [get_file(time) for time in pope_times]

# Calculate PDF
mean_w = np.ma.zeros(num_levels)
median_w = np.ma.zeros(num_levels)
ninety_w = np.ma.zeros(num_levels)
ninety_five_w = np.ma.zeros(num_levels)
ninety_nine_w = np.ma.zeros(num_levels)
mean_z = np.ma.zeros(num_levels)
median_z = np.ma.zeros(num_levels)
ninety_z = np.ma.zeros(num_levels)
ninety_five_z = np.ma.zeros(num_levels)
ninety_nine_z = np.ma.zeros(num_levels)
bins = np.arange(-10,40,1)
bins_z = np.arange(0,60,1)
w_hist = np.ma.zeros((num_levels, len(bins)-1))
w_cfad = np.ma.zeros((num_levels, len(bins)-1))
z_hist = np.ma.zeros((num_levels, len(bins_z)-1))
ws = [da.from_delayed(arrays, shape=(num_levels,81,111),
                      dtype=float) for arrays in ws]
ws = da.stack(ws, axis=0)
 
import time
for levels in range(0,num_levels):
    t1 = time.time()  
    w_level = ws[:,levels,:,:]   
    print(str(levels) + '/' + str(num_levels))
    array_shape = w_level.shape
    # Take chunks of w_level and remove NaNs
    w_new = []
    for i in range(0, array_shape[0], int(array_shape[0]/10)):
        print('Processing chunk ' + str(i))
        w_chunk = np.array(w_level[i:(i+int(array_shape[0]/10)),:,:])
        w_chunk = w_chunk[~np.isnan(w_chunk)]
        w_new.append(w_chunk.flatten())
    
    w_new = np.concatenate(w_new)
    w_new = w_new.flatten()
    print(w_new)
    num_elems = array_shape[0]*array_shape[1]*array_shape[2]    
    print(w_new.shape)
    print(mean_w.shape)
    if(len(w_new) > 0):
        mean_w[levels] = np.nanmean(w_new)  
        median_w[levels] = np.nanpercentile(w_new, 50)
        ninety_w[levels] = np.nanpercentile(w_new, 90)
        ninety_five_w[levels] = np.percentile(w_new, 95)
        ninety_nine_w[levels] = np.percentile(w_new, 99)
    else:
        mean_w[levels] = np.nan
        median_w[levels] = np.nan
        ninety_w[levels] = np.nan
        ninety_five_w[levels] = np.nan
        ninety_nine_w[levels] = np.nan
    t2 = time.time() - t1
    print('Total time in s: ' + str(t2))
    print('Time per scan = ' + str(t2/array_shape[0]))
        
print('Writing netCDF file...')
# Save to netCDF file
out_netcdf = Dataset('wpdfregime' + str(pope_regime) + '.cdf', 'w')
out_netcdf.createDimension('levels', num_levels)
mean_file = out_netcdf.createVariable('mean', mean_w.dtype, ('levels',))
mean_file.long_name = 'Mean w'
mean_file.units = 'm s-1'
mean_file[:] = mean_w

median_file = out_netcdf.createVariable('median', median_w.dtype, ('levels',))
median_file.long_name = 'median w'
median_file.units = 'm s-1'
median_file[:] = median_w

ninety_file = out_netcdf.createVariable('ninety', ninety_w.dtype, ('levels',))
ninety_file.long_name = '90% w'
ninety_file.units = 'm s-1'
ninety_file[:] = ninety_w

n5_file = out_netcdf.createVariable('ninety_five', ninety_five_w.dtype, ('levels',))
n5_file.long_name = '95W w'
n5_file.units = 'm s-1'
n5_file[:] = ninety_five_w

n5_file = out_netcdf.createVariable('ninety_nine', ninety_five_w.dtype, ('levels',))
n5_file.long_name = '99W w'
n5_file.units = 'm s-1'
n5_file[:] = ninety_nine_w

z_file = out_netcdf.createVariable('z', ninety_five_w.dtype, ('levels',))
z_file.long_name = 'z'
z_file.units = 'm'
z_file[:] = z_levels
                                       
out_netcdf.close()
