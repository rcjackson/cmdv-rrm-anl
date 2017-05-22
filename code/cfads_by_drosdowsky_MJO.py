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
import pandas

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
mjo_index_file = '/home/rjackson/data/rmm.74toRealtime.txt'

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

def get_z(time):
    pyart_grid = get_grid_from_dda(time)
    z = pyart_grid.fields['reflectivity']['data']
    z_individual = []
    level_individual = []

    # Find deep convective cores and get max updraft speeds
    for levels in range(0,num_levels-1): 
        ref = z[levels]
        ref = ref[ref > 1]
        z_individual.append(ref.flatten())
        level_individual.append(levels*np.ones(len(ref.flatten())))
                        
    # Convert to list of individual max w's for each updraft
    if(len(z_individual) > 0):
        z_individual = np.concatenate(z_individual)
        level_individual = np.concatenate(level_individual)

    return_array = np.ma.zeros((len(z_individual),2))
    return_array[:,0] = z_individual
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

# Load MJO index
data = pandas.read_csv(mjo_index_file,
                       header=2,
                       delim_whitespace=True)
data_matrix = np.ma.array(data.values)
year = data_matrix[:,0]
month = data_matrix[:,1]
day = data_matrix[:,2]
index = data_matrix[:,5]
mjodates = []
for i in range(0, len(day)):
    mjodates.append(datetime(year=int(year[i]),
                             month=int(month[i]),
                             day=int(day[i])))

# Since grids are uniform, calculate beam crossing angle for first grid and
# apply to all
first_grid = get_grid_from_dda(times[0])
bca = get_bca(first_grid) 
num_levels = 40
z_levels = np.arange(0.5,0.5*(num_levels+1),0.5)*1000
count = 0
dros_regime = int(sys.argv[2])
mjo_index_start = int(sys.argv[3])
mjo_index_end = int(sys.argv[4])

# Filter out data not in Pope regime
dros_times = []
for time in times:
    # Look for date in Pope regime data
    cur_date = datetime(year=time.year, month=time.month, day=time.day)
    inds = np.where([day <= cur_date for day in drosdates])
    inds_mjo = np.where([day <= cur_date for day in mjodates])
    dros_index = inds[0][-1]
    inds_mjo = inds_mjo[0][-1]
    if(groups[dros_index] == dros_regime and 
       index[inds_mjo] >= mjo_index_start and 
       index[inds_mjo] <= mjo_index_end):
        print((drosdates[dros_index], time))
        dros_times.append(time)

in_netcdf.close()

# Get delayed structure to load files in parallel
get_file = delayed(get_z)

# Calculate CFAD
bin_spacing = 1
bins_z = np.arange(0,60,bin_spacing)
cfad = np.ma.zeros((num_levels, len(bins_z)-1))
print('Doing parallel grid loading...')
import time
t1 = time.time()
ws = []

for i in range(0, len(dros_times), int(len(dros_times)/8)):
    ws_temp = [get_file(times) for times in dros_times[i:(i+len(dros_times)/8)]]
    ws_temp = compute(*ws_temp)
    # Do histogram for each chunk
    
    level_individual = ws_temp[0][:,1] 
    z_individual = ws_temp[0][:,0]
    for levels in range(0,num_levels):
        z_new = z_individual[level_individual == levels] 
        if(len(z_new) > 0):
            hist, bins = np.histogram(z_new, bins=bins_z, normed=False)
            cfad[levels,:] = cfad[levels,:] + hist

# Normalize final histogram
for levels in range(0,num_levels):
    cfad[levels,:] = cfad[levels,:]/np.sum(cfad[levels,:])/bin_spacing

t2 = time.time() - t1
print('Total time in s: ' + str(t2))
print('Time per scan = ' + str(t2/len(dros_times)))

# Save to netCDF file
out_netcdf = Dataset(('cfadregime' + 
                       str(dros_regime) + 
                       '_dros' + 
                       str(mjo_index_start) +
                       str(mjo_index_end) + '.cdf'), 'w')
out_netcdf.createDimension('levels', num_levels)
out_netcdf.createDimension('bins', len(bins_z)-1)

cfad_file = out_netcdf.createVariable('cfad', cfad.dtype, ('levels', 'bins'))
cfad_file.long_name = 'CFAD of Reflectivity'
cfad_file.units = 'normalized'
cfad_file[:] = cfad

z_file = out_netcdf.createVariable('z', z_levels.dtype, ('levels',))
z_file.long_name = 'z'
z_file.units = 'm'
z_file[:] = z_levels

bins_z_file = out_netcdf.createVariable('bins_z', bins_z.dtype, ('bins',))
bins_z_file.long_name = 'z'
bins_z_file.units = 'm'
bins_z_file[:] = (bins_z[1:]+bins_z[0:-1])/2
                                       
out_netcdf.close()
