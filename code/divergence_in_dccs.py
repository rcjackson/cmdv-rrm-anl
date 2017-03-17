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
from dask import delayed, compute
import time
import sys
from scipy import ndimage

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
    z = pyart_grid.fields['reflectivity']['data']
    bca = np.ma.masked_invalid(bca)
    u = pyart_grid.fields['eastward_wind']['data']
    v = pyart_grid.fields['northward_wind']['data']   
    u[u.mask == True] = np.nan
    v[v.mask == True] = np.nan                            
    divergence = np.gradient(u, 1e6, axis=2) + np.gradient(v, 1e6, axis=1)
    for levels in range(0,num_levels-1):
        w_outside_updraft = np.logical_or(w[levels] < 1, w[levels] > 99.0)
        outside_dd_lobes = np.logical_or(bca < math.pi/6, bca > 5*math.pi/6)
        w[levels] = np.ma.masked_where(np.logical_or(w_outside_updraft,
                                                     outside_dd_lobes), w[levels])
        z[levels] = np.ma.masked_where(np.logical_or(w_outside_updraft,
                                                     outside_dd_lobes), z[levels])
        divergence[levels] = np.ma.masked_where(np.logical_or(w_outside_updraft,
                                                              outside_dd_lobes), 
                                                divergence[levels]) 

    grid_z = pyart_grid.point_z['data']

    # Set mask to exclude data outside of updrafts
    w_temp = deepcopy(w)
    w_temp[~w_temp.mask] = 1
    w_temp[w_temp.mask] = 0
    w_temp.mask = False

    six_connected_structure = [[[0,0,0],
                                [0,1,0],
                                [0,0,0]],
                               [[0,1,0],
                                [1,1,1],
                                [0,1,0]],
                               [[0,0,0],
                                [0,1,0],
                                [0,0,0]]]
    updrafts, num_updrafts = ndimage.measurements.label(w_temp, 
                                                        structure=six_connected_structure)
    
    # Get statistics in continous regions
    index = np.arange(0, num_updrafts + 1)
    max_z = ndimage.measurements.maximum(grid_z, 
                                         labels=updrafts, 
                                         index=index)
    min_z = ndimage.measurements.minimum(grid_z, 
                                         labels=updrafts,
                                         index=index)
    
    div_individual = []
    level_individual = []
    label_individual = []
    count_individual = []

    # Find deep convective cores and get max updraft speeds
    for levels in range(0,num_levels-1):
        label_level = updrafts[levels]
        masked_array = np.ma.zeros(updrafts.shape)
        masked_array.mask = True
        w_temp = w[levels]
        div = divergence[levels]
        for labels in range(1, len(max_z)-1):
            indicies = np.ma.where(label_level == labels)                                
                        
            if(len(indicies[0]) > 0  
               and max_z[labels] >= 6000
               and min_z[labels] <= 1000):
                div_individual.append(div[indicies])
                level_individual.append(levels*np.ones(len(div[indicies])))
                label_individual.append(labels*np.ones(len(div[indicies])))
                count_individual.append(count*np.ones(len(div[indicies])))
                
    # Convert to list of individual max w's for each updraft
    if(len(div_individual) > 0):
        div_individual = np.concatenate(div_individual)
        level_individual = np.concatenate(level_individual)

    # Very large vertical velocities aloft
    #if(len(max_w_individual) > 0):
    #    if(np.max(max_w_individual) > 60):
    #        print('Very large vertical velocity found:')
    #        print(time)
    #        max_w_individual = np.array([])
    #        level_individual = np.array([])
    print(div_individual)
    print(level_individual)
    return_array = np.ma.zeros((len(div_individual),2))
    return_array[:,0] = div_individual
    return_array[:,1] = level_individual
    return return_array
       

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
    if(groups[pope_index] == pope_regime):
        print((popedates[pope_index], time))
        pope_times.append(time)

in_netcdf.close()

# Get delayed structure to load files in parallel
get_file = delayed(get_updrafts)

# Calculate PDF
mean_w = np.ma.zeros(num_levels)
median_w = np.ma.zeros(num_levels)
ninety_w = np.ma.zeros(num_levels)
ninety_five_w = np.ma.zeros(num_levels)
ninety_nine_w = np.ma.zeros(num_levels)
five_w = np.ma.zeros(num_levels)
one_w = np.ma.zeros(num_levels)
ten_w = np.ma.zeros(num_levels)

bins = np.arange(-10,40,1)
bins_z = np.arange(0,60,1)
print('Doing parallel grid loading...')
import time
t1 = time.time()
ws = []
for i in range(0, len(pope_times), int(len(pope_times)/4)):
    ws_temp = [get_file(times) for times in pope_times[i:(i+len(pope_times)/4)]]
    ws_temp = compute(*ws_temp)
    ws.append(ws_temp)

for arrays in ws:
    array_temp = np.concatenate(arrays)
    print(array_temp.shape) 

ws = np.concatenate([np.concatenate(arrays) for arrays in ws])

t2 = time.time() - t1
print('Total time in s: ' + str(t2))
print('Time per scan = ' + str(t2/len(pope_times)))
level_individual = ws[:,1] 
w_individual = ws[:,0]
print(len(level_individual))
print(len(w_individual))
for levels in range(0,num_levels):      
    w_new = w_individual[level_individual == levels]
    if(len(w_new) > 0):
        mean_w[levels] = np.nanmean(w_new)  
        median_w[levels] = np.nanpercentile(w_new, 50)
        ninety_w[levels] = np.nanpercentile(w_new, 90)
        ninety_five_w[levels] = np.nanpercentile(w_new, 95)
        ninety_nine_w[levels] = np.nanpercentile(w_new, 99)
        five_w[levels] = np.nanpercentile(w_new, 5)
        one_w[levels] = np.nanpercentile(w_new, 1)
        ten_w[levels] = np.nanpercentile(w_new, 10)
    else:
        mean_w[levels] = np.nan
        median_w[levels] = np.nan
        ninety_w[levels] = np.nan
        ninety_five_w[levels] = np.nan
        ninety_nine_w[levels] = np.nan
        five_w[levels] = np.nan
        one_w[levels] = np.nan
        ten_w[levels] = np.nan
    
print('Writing netCDF file...')
# Save to netCDF file
out_netcdf = Dataset('divpdfregime' + str(pope_regime) + '_varble.cdf', 'w')
out_netcdf.createDimension('levels', num_levels)
mean_file = out_netcdf.createVariable('mean', mean_w.dtype, ('levels',))
mean_file.long_name = 'Mean div'
mean_file.units = 'm s-1'
mean_file[:] = mean_w

median_file = out_netcdf.createVariable('median', median_w.dtype, ('levels',))
median_file.long_name = 'median div'
median_file.units = 'm s-1'
median_file[:] = median_w

ninety_file = out_netcdf.createVariable('ninety', ninety_w.dtype, ('levels',))
ninety_file.long_name = '90% div'
ninety_file.units = 'm s-1'
ninety_file[:] = ninety_w

n5_file = out_netcdf.createVariable('ninety_five', ninety_five_w.dtype, ('levels',))
n5_file.long_name = '95% div'
n5_file.units = 'm s-1'
n5_file[:] = ninety_five_w

n5_file = out_netcdf.createVariable('ninety_nine', ninety_five_w.dtype, ('levels',))
n5_file.long_name = '99% div'
n5_file.units = 'm s-1'
n5_file[:] = ninety_nine_w

one_file = out_netcdf.createVariable('one', ninety_w.dtype, ('levels',))
one_file.long_name = '1% div'
one_file.units = 'm s-1'
one_file[:] = one_w

five_file = out_netcdf.createVariable('five', ninety_five_w.dtype, ('levels',))
five_file.long_name = '5% div'
five_file.units = 'm s-1'
five_file[:] = five_w

ten_file = out_netcdf.createVariable('ten', ninety_five_w.dtype, ('levels',))
ten_file.long_name = '10% div'
ten_file.units = 'm s-1'
ten_file[:] = ten_w

z_file = out_netcdf.createVariable('z', ninety_five_w.dtype, ('levels',))
z_file.long_name = 'z'
z_file.units = 'm'
z_file[:] = z_levels
                                       
out_netcdf.close()
