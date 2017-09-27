import matplotlib
matplotlib.use('Agg')
import pyart
import numpy as np
import glob
import math
import dask.bag as db
import time
import sys

from datetime import datetime, timedelta
from copy import deepcopy
from distributed import Client, LocalCluster
from dask import delayed, compute
from netCDF4 import Dataset
from scipy import ndimage

# Start a cluster with x workers
client = Client(scheduler_file=sys.argv[1])
land_ocean = int(sys.argv[3])

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
    
    # Parse all of the dates and time in the 
    # interval and add them to the time list
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
        
    
    # Sort time list and make sure time are 
    # at least xx min apart
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

def get_updrafts(time, bca, ocean_mask, land_ocean):
    try:
        pyart_grid = get_grid_from_dda(time)
        w = pyart_grid.fields['upward_air_velocity']['data']
    except:
        return_array = np.ma.zeros((0,3))
        max_w_individual = np.array([])
        level_individual = np.array([])
        flux_individual = np.array([])
        return_array[:,0] = max_w_individual
        return_array[:,1] = level_individual
        return_array[:,2] = flux_individual
        return return_array

    z = pyart_grid.fields['reflectivity']['data']
    bca = np.ma.masked_invalid(bca)
    
    if(land_ocean == 1):
        truth_value = True
    else:
        truth_value = False

    for levels in range(0, num_levels-1):
        w_outside_updraft = np.logical_or(w[levels] < 1, w[levels] > 99.0)
        outside_dd_lobes = np.logical_or(bca < math.pi/6, bca > 5*math.pi/6)
        if(land_ocean < 2):
            w[levels] = np.ma.masked_where(np.logical_or(ocean_mask[levels] == truth_value,
                np.logical_or(w_outside_updraft, outside_dd_lobes)), w[levels])
            z[levels] = np.ma.masked_where(np.logical_or(ocean_mask[levels] == truth_value,
                np.logical_or(w_outside_updraft, outside_dd_lobes)), w[levels])
        else:
            w[levels] = np.ma.masked_where(np.logical_or(
                w_outside_updraft, outside_dd_lobes), w[levels])
            z[levels] = np.ma.masked_where(np.logical_or(
                w_outside_updraft, outside_dd_lobes), w[levels])
       
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
    flux_individual = []

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
                mflux = np.ma.sum(w_temp[indicies])*1.0*np.exp(levels*0.5/8)*1e6
                flux_individual.append(mflux) 
                level_individual.append(levels)
                label_individual.append(labels)
                count_individual.append(count)
    
    # Convert to list of individual max w's for each updraft
    max_w_individual = np.array(max_w_individual)
    level_individual = np.array(level_individual)
    flux_individual = np.array(flux_individual)

    # Very large vertical velocities aloft
    if(len(max_w_individual) > 0):
        if(np.max(max_w_individual) > 60):
            print('Very large vertical velocity found:')
            print(time)
            max_w_individual = np.array([])
            level_individual = np.array([])
            flux_individual = np.array([])
    return_array = np.ma.zeros((len(max_w_individual),3))
    return_array[:,0] = max_w_individual
    return_array[:,1] = level_individual
    return_array[:,2] = flux_individual
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
first_grid = pyart.io.read_grid('ocean_mask.nc')
bca = get_bca(first_grid) 
num_levels = 40
z_levels = np.arange(0.5,0.5*(num_levels+1),0.5)*1000
count = 0
dros_regime = int(sys.argv[2])
ocean_mask = first_grid.fields['ocean_mask']['data']
ocean_mask = ocean_mask.mask

# Filter out data not in Drosdowsky regime
dros_times = []
for time in times:
    # Look for date in Drosdowsky regime data
    cur_date = datetime(year=time.year, month=time.month, day=time.day)
    inds = np.where([day <= cur_date for day in drosdates])
    dros_index = inds[0][-1]
    if(groups[dros_index] == dros_regime):
        dros_times.append(time)

in_netcdf.close()

get_updrafts_func = lambda x: get_updrafts(
    x, bca, ocean_mask, land_ocean)

# Get delayed structure to load files in parallel
#get_file = delayed(get_updrafts_func)

# Calculate PDF
mean_w = np.ma.zeros(num_levels)
median_w = np.ma.zeros(num_levels)
ninety_w = np.ma.zeros(num_levels)
ninety_five_w = np.ma.zeros(num_levels)
ninety_nine_w = np.ma.zeros(num_levels)
mean_f = np.ma.zeros(num_levels)
median_f = np.ma.zeros(num_levels)
ninety_f = np.ma.zeros(num_levels)
ninety_five_f = np.ma.zeros(num_levels)
ninety_nine_f = np.ma.zeros(num_levels)
bins = np.arange(-10,40,1)
bins_z = np.arange(0,60,1)
print('Doing parallel grid loading...')
import time
t1 = time.time()
ws = []
print(len(dros_times))
the_bag = db.from_sequence(dros_times)
ws = the_bag.map(get_updrafts_func).compute()

print(ws)
ws = np.concatenate(ws, axis=0)
t2 = time.time() - t1
print('Total time in s: ' + str(t2))
print('Time per scan = ' + str(t2/len(dros_times)))
level_individual = ws[:,1] 
w_individual = ws[:,0]
flux_individual = ws[:,2]
for levels in range(0,num_levels):      
    w_new = w_individual[level_individual == levels]
    flux = flux_individual[level_individual == levels]
    if(len(w_new) > 0):
        mean_w[levels] = np.nanmean(w_new)  
        median_w[levels] = np.nanpercentile(w_new, 50)
        ninety_w[levels] = np.nanpercentile(w_new, 90)
        ninety_five_w[levels] = np.percentile(w_new, 95)
        ninety_nine_w[levels] = np.percentile(w_new, 99)
        mean_f[levels] = np.nanmean(flux)  
        median_f[levels] = np.nanpercentile(flux, 50)
        ninety_f[levels] = np.nanpercentile(flux, 90)
        ninety_five_f[levels] = np.percentile(flux, 95)
        ninety_nine_f[levels] = np.percentile(flux, 99)
    else:
        mean_w[levels] = np.nan
        median_w[levels] = np.nan
        ninety_w[levels] = np.nan
        ninety_five_w[levels] = np.nan
        ninety_nine_w[levels] = np.nan
        mean_f[levels] = np.nan
        median_f[levels] = np.nan
        ninety_f[levels] = np.nan
        ninety_five_f[levels] = np.nan
        ninety_nine_f[levels] = np.nan

print('Writing netCDF file...')

# Save to netCDF file
out_netcdf = Dataset(('wpdfregime_dros' + str(dros_regime) + 
                      '_varble' + str(land_ocean) + '.cdf'), 'w')
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

mean_file = out_netcdf.createVariable('mean_f', mean_w.dtype, ('levels',))
mean_file.long_name = 'mean mass flux'
mean_file.units = 'm s-1'
mean_file[:] = mean_f

median_file = out_netcdf.createVariable('median_f', median_w.dtype, ('levels',))
median_file.long_name = 'median mass flux'
median_file.units = 'kg s-1'
median_file[:] = median_f

ninety_file = out_netcdf.createVariable('ninety_f', ninety_w.dtype, ('levels',))
ninety_file.long_name = '90% mass flux'
ninety_file.units = 'kg s-1'
ninety_file[:] = ninety_f

n5_file = out_netcdf.createVariable('ninety_five_f', ninety_five_w.dtype, ('levels',))
n5_file.long_name = '95% mass flux'
n5_file.units = 'kg s-1'
n5_file[:] = ninety_five_f

n5_file = out_netcdf.createVariable('ninety_nine_f', ninety_five_w.dtype, ('levels',))
n5_file.long_name = '99 mass flux'
n5_file.units = 'kg s-1'
n5_file[:] = ninety_nine_f
                           
out_netcdf.close()
