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
start_year = 2006
start_day = 1
start_month = 1
start_hour = 1
start_minute = 0
start_second = 0

end_year = 2006
end_month = 3
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

    for levels in range(0,num_levels-1):
        w_outside_updraft = np.logical_or(w[levels] < 1, w[levels] > 99.0)
        outside_dd_lobes = np.logical_or(bca < math.pi/6, bca > 5*math.pi/6)
        w[levels] = np.ma.masked_where(np.logical_or(w_outside_updraft,
                                                     outside_dd_lobes), w[levels])
        z[levels] = np.ma.masked_where(np.logical_or(w_outside_updraft,
                                                     outside_dd_lobes), z[levels])
       
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
    index=np.arange(0, num_updrafts + 1)
    max_z = ndimage.measurements.maximum(grid_z, 
                                         labels=updrafts, 
                                         index=index)
    min_z = ndimage.measurements.minimum(grid_z, 
                                         labels=updrafts,
                                         index=index)
    
    max_w_individual = []
    level_individual = []
    label_individual = []
    count_individual = []

    # Find deep convective cores and get max updraft speeds
    for levels in range(0,num_levels-1):
        label_level = updrafts[levels]
        masked_array = np.ma.zeros(updrafts.shape)
        masked_array.mask = True
        w_temp = w[levels]
        
        for labels in range(1, len(max_z)-1):
            indicies = np.ma.where(label_level == labels)                                
                        
            if(len(indicies[0]) > 0  
               and max_z[labels] >= 15000
               and min_z[labels] <= 1000):
                max_w_individual.append(max(w_temp[indicies]))
                level_individual.append(levels)
                label_individual.append(labels)
                count_individual.append(count)
    
    # Convert to list of individual max w's for each updraft
    max_w_individual = np.array(max_w_individual)
    level_individual = np.array(level_individual)

    # Very large vertical velocities aloft
    if(len(max_w_individual) > 0):
        if(np.max(max_w_individual) > 60):
            print('Very large vertical velocity found:')
            print(time)
            max_w_individual = np.array([])
            level_individual = np.array([])
    return_array = np.ma.zeros((len(max_w_individual),2))
    return_array[:,0] = max_w_individual
    return_array[:,1] = level_individual
    return return_array
    
# Plot the radars from given time.
times = get_dda_times(start_year, start_month, start_day,
                      start_hour, start_minute, end_year,
                      end_month, end_day, end_hour, 
                      end_minute, minute_interval=0)

in_netcdf = Dataset('/home/rjackson/data/Drosdowsky.cdf', 
                    mode='r')            
year = in_netcdf.variables['year'][:]
month = in_netcdf.variables['month'][:]
day = in_netcdf.variables['day'][:]
groups = in_netcdf.variables['groups'][:]

drosdates = []
for i in range(0,len(day)):
    drosdates.append(datetime(year=int(year[i]),
                              month=int(month[i]),
                              day=int(day[i])))

# Since grids are uniform, calculate beam crossing angle for first grid and
# apply to all
num_levels = 40
z_levels = np.arange(0.5,0.5*(num_levels+1),0.5)*1000
count = 0
dros_regime = int(sys.argv[2])

# Filter out data not in Drosdowsky regime
dros_times = []
for time in times:
    # Look for date in Drosdwosky regime data
    cur_date = datetime(year=time.year, month=time.month, day=time.day)
    inds = np.where([day <= cur_date for day in drosdates])
    dros_index = inds[0][-1]
    if(groups[dros_index] == dros_regime):
        print((drosdates[dros_index], time))
        dros_times.append(time)

in_netcdf.close()

# Load SCAI
CY_cdf = Dataset('/home/rjackson/data/num_clusters.cdf', mode='r')
classification = CY_cdf.variables['N'][:]
d1 = CY_cdf.variables['d1'][:]
year_t = CY_cdf.variables['year'][:]
month_t = CY_cdf.variables['month'][:]
day_t = CY_cdf.variables['day'][:]
hour_t = CY_cdf.variables['hour'][:]
minute_t = CY_cdf.variables['minute'][:]

L = 350.0/5
a = 5.0
Nmax = pow((L/a),2)
SCAI = ((classification))/Nmax*d1/(L)*1000

# Get times from Tobin et al. (2012) file
date_array_tobin = []
for i in range(0,len(year_t)):
    dat = datetime(year_t[i],month_t[i],day_t[i],hour_t[i],minute_t[i],)
    date_array_tobin.append(dat)

# Sort times into organized (SCAI < 40) and disorganized (SCAI > 40)
dros_times_disorganized = []
dros_times_organized = []
for times in dros_times:
    # Look for date in Drosdwosky regime data
    cur_date = datetime(year=times.year, 
                        month=times.month, 
                        day=times.day,
                        hour=times.hour)
    inds = np.where([day <= cur_date for day in date_array_tobin])
    tobin_index = inds[0][-1]
    if(SCAI[tobin_index] > 40):
        dros_times_disorganized.append(times)
    elif(SCAI[tobin_index] <= 40):
        dros_times_organized.append(times)
 
# Get delayed structure to load files in parallel
get_file = delayed(get_updrafts)

# Calculate PDF
mean_w_organized = np.ma.zeros(num_levels)
median_w_organized = np.ma.zeros(num_levels)
ninety_w_organized = np.ma.zeros(num_levels)
ninety_five_w_organized = np.ma.zeros(num_levels)
ninety_nine_w_organized = np.ma.zeros(num_levels)
mean_w_disorganized = np.ma.zeros(num_levels)
median_w_disorganized = np.ma.zeros(num_levels)
ninety_w_disorganized = np.ma.zeros(num_levels)
ninety_five_w_disorganized = np.ma.zeros(num_levels)
ninety_nine_w_disorganized = np.ma.zeros(num_levels)
bins = np.arange(-10,40,1)
bins_z = np.arange(0,60,1)
print('Doing parallel grid loading...')

import time

# Group updrafts in organized convection together
t1 = time.time()
ws_organized = []
num_times = len(dros_times_organized)
for i in range(0, num_times, int(num_times/4)):
    ws_temp = [get_file(times) 
               for times in dros_times_organized[i:(i+num_times/4)]]
    ws_temp = compute(*ws_temp)
    ws_organized.append(ws_temp)
 
# Group updrafts in disorganized convection together
t1 = time.time()
num_times = len(dros_times_disorganized)
ws_disorganized = []
for i in range(0, num_times, int(num_times/4)):
    ws_temp = [get_file(times) 
               for times in dros_times_disorganized[i:(i+num_times/4)]]
    ws_temp = compute(*ws_temp)
    ws_disorganized.append(ws_temp)

ws_organized = np.concatenate([np.concatenate(arrays) 
                               for arrays in ws_organized])
ws_disorganized = np.concatenate([np.concatenate(arrays) 
                                  for arrays in ws_disorganized])
t2 = time.time() - t1
print('Total time in s: ' + str(t2))
print('Time per scan = ' + str(t2/len(dros_times)))

level_individual_organized = ws_organized[:,1] 
w_individual_organized = ws_organized[:,0]
level_individual_disorganized = ws_disorganized[:,1] 
w_individual_disorganized = ws_disorganized[:,0]
for levels in range(0,num_levels):      
    w_new = w_individual_organized[level_individual_organized == levels]
    if(len(w_new) > 0):
        mean_w_organized[levels] = np.nanmean(w_new)  
        median_w_organized[levels] = np.nanpercentile(w_new, 50)
        ninety_w_organized[levels] = np.nanpercentile(w_new, 90)
        ninety_five_w_organized[levels] = np.percentile(w_new, 95)
        ninety_nine_w_organized[levels] = np.percentile(w_new, 99)
    else:
        mean_w_organized[levels] = np.nan
        median_w_organized[levels] = np.nan
        ninety_w_organized[levels] = np.nan
        ninety_five_w_organized[levels] = np.nan
        ninety_nine_w_organized[levels] = np.nan
    
    w_new = w_individual_disorganized[level_individual_disorganized == levels]
    if(len(w_new) > 0):
        mean_w_disorganized[levels] = np.nanmean(w_new)  
        median_w_disorganized[levels] = np.nanpercentile(w_new, 50)
        ninety_w_disorganized[levels] = np.nanpercentile(w_new, 90)
        ninety_five_w_disorganized[levels] = np.percentile(w_new, 95)
        ninety_nine_w_disorganized[levels] = np.percentile(w_new, 99)
    else:
        mean_w_disorganized[levels] = np.nan
        median_w_disorganized[levels] = np.nan
        ninety_w_disorganized[levels] = np.nan
        ninety_five_w_disorganized[levels] = np.nan
        ninety_nine_w_disorganized[levels] = np.nan
                    
print('Writing netCDF file...')

# Save to netCDF file
out_netcdf = Dataset('wpdfregime_dros' + str(dros_regime) + '_varble_org.cdf', 'w')
out_netcdf.createDimension('levels', num_levels)
mean_file = out_netcdf.createVariable('mean_organized', 
                                      mean_w_organized.dtype, 
                                      ('levels',))
mean_file.long_name = 'Mean w'
mean_file.units = 'm s-1'
mean_file[:] = mean_w_organized

median_file = out_netcdf.createVariable('median_organized', 
                                        median_w_organized.dtype, 
                                        ('levels',))
median_file.long_name = 'median w'
median_file.units = 'm s-1'
median_file[:] = median_w_organized

ninety_file = out_netcdf.createVariable('ninety_organized', 
                                        ninety_w_organized.dtype, 
                                        ('levels',))
ninety_file.long_name = '90% w'
ninety_file.units = 'm s-1'
ninety_file[:] = ninety_w_organized

n5_file = out_netcdf.createVariable('ninety_five_organized', 
                                     ninety_five_w_organized.dtype, 
                                     ('levels',))
n5_file.long_name = '95W w'
n5_file.units = 'm s-1'
n5_file[:] = ninety_five_w_organized

n5_file = out_netcdf.createVariable('ninety_nine_organized',
                                    ninety_five_w_organized.dtype, 
                                    ('levels',))
n5_file.long_name = '99W w'
n5_file.units = 'm s-1'
n5_file[:] = ninety_nine_w_organized

mean_file = out_netcdf.createVariable('mean_disorganized', 
                                      mean_w_disorganized.dtype, 
                                      ('levels',))
mean_file.long_name = 'Mean w'
mean_file.units = 'm s-1'
mean_file[:] = mean_w_disorganized

median_file = out_netcdf.createVariable('median_disorganized', 
                                         median_w_disorganized.dtype, 
                                        ('levels',))
median_file.long_name = 'median w'
median_file.units = 'm s-1'
median_file[:] = median_w_disorganized

ninety_file = out_netcdf.createVariable('ninety_disorganized', 
                                        ninety_w_disorganized.dtype, 
                                        ('levels',))
ninety_file.long_name = '90% w'
ninety_file.units = 'm s-1'
ninety_file[:] = ninety_w_disorganized

n5_file = out_netcdf.createVariable('ninety_five_disorganized', 
                                    ninety_five_w_disorganized.dtype,
                                    ('levels',))
n5_file.long_name = '95W w'
n5_file.units = 'm s-1'
n5_file[:] = ninety_five_w_disorganized

n5_file = out_netcdf.createVariable('ninety_nine_disorganized', 
                                    ninety_five_w_disorganized.dtype, 
                                    ('levels',))
n5_file.long_name = '99W w'
n5_file.units = 'm s-1'
n5_file[:] = ninety_nine_w_disorganized

z_file = out_netcdf.createVariable('z', 
                                   ninety_five_w_organized.dtype, 
                                   ('levels',))
z_file.long_name = 'z'
z_file.units = 'm'
z_file[:] = z_levels
                                       
out_netcdf.close()
