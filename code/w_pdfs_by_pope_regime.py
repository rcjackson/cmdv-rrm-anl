import pyart
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
from copy import deepcopy
import glob
import math

# Input the range of dates and time wanted for the collection of images
start_year = 2005
start_day = 1
start_month = 1
start_hour = 1
start_minute = 0
start_second = 0

end_year = 2011
end_month = 5
end_day = 24
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

# Plot the radars from given time.

times = get_dda_times(start_year, start_month, start_day,
                      start_hour, start_minute, end_year,
                      end_month, end_day, end_hour, 
                      end_minute, minute_interval=0)

in_netcdf = Dataset('/lcrc/group/earthscience/rjackson/data/Pope_regime.cdf', mode='r')            
year = in_netcdf.variables['year'][:]
month = in_netcdf.variables['month'][:]
day = in_netcdf.variables['day'][:]
groups = in_netcdf.variables['groups'][:]

popedates = []
for i in range(0,len(day)):
    popedates.append(datetime(year=int(year[i]),
                              month=int(month[i]),
                              day=int(day[i])))

num_levels = 40
z_levels = np.arange(0.5,0.5*(num_levels+1),0.5)*1000
count = 0
pope_regime = 0
dc = np.ma.zeros((len(times), 81, 111))
ws = np.ma.zeros((len(times), num_levels, 81, 111))
ws_all = np.ma.zeros((len(times), num_levels, 81, 111))
rfs = np.ma.zeros((len(times), num_levels, 81, 111))

for time in times:
    # Look for date in Pope regime data
    cur_date = datetime(year=time.year, month=time.month, day=time.day)
    inds = np.where([day <= cur_date for day in popedates])
    pope_index = inds[0][-1]
    print((popedates[pope_index], time))
    if(groups[pope_index] == pope_regime):
        pyart_grid = get_grid_from_dda(time)
        bca = get_bca(pyart_grid)
        w = pyart_grid.fields['upward_air_velocity']['data']
        updraft_depth = np.zeros(w[0].shape)
        z = pyart_grid.fields['reflectivity']['data']
        for levels in range(0,num_levels-1):
            w[levels] = np.ma.masked_where(np.logical_or(np.logical_or(bca < math.pi/6,
                                                                       bca > 5*math.pi/6), 
                                                         z[levels] < 1), w[levels])

            is_in_updraft = w[levels] > 1
            is_in_updraft_next = w[levels+1] > 1
            both_in_updraft = np.logical_or(np.logical_and(is_in_updraft,
                                                           is_in_updraft_next),
                                            updraft_depth > 10)
        

            add_one = np.where(both_in_updraft)
            set_to_zero = np.where(~both_in_updraft)
            if(len(add_one[0]) > 0):
                updraft_depth[add_one[0], add_one[1]] = updraft_depth[add_one[0], 
                                                                      add_one[1]] + 1 
                updraft_depth[set_to_zero[0], set_to_zero[1]] = 0       
    
        dc[count] = updraft_depth
        ws[count] = w
        ws_all[count] = w
        rfs[count] = z
        count = count + 1

in_netcdf.close()

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
dims = ws.shape
bins = np.arange(-10,40,1)
bins_z = np.arange(0,60,1)
w_hist = np.ma.zeros((num_levels, len(bins)-1))
w_cfad = np.ma.zeros((num_levels, len(bins)-1))
z_hist = np.ma.zeros((num_levels, len(bins_z)-1))

for levels in range(0,num_levels):
    w_level = ws[:,levels,:,:]
    r_level = rfs[:,levels,:,:]
    for i in range(0, dims[0]):
        w_level[i,:,:] = np.ma.masked_where(np.logical_or(dc[i,:,:] < 10,
                                                          r_level[i,:,:] < 1),
                                                          w_level[i,:,:])
        
    ws_in_core = w_level[~w_level.mask]
    zs_in_core = r_level[~w_level.mask]
    mean_w[levels] = np.ma.mean(ws_in_core[ws_in_core < 99])
    median_w[levels] = np.ma.median(ws_in_core[ws_in_core < 99])
    mean_z[levels] = np.ma.mean(zs_in_core[ws_in_core < 99])
    median_z[levels] = np.ma.median(zs_in_core[ws_in_core < 99])
    counts, bins = np.histogram(ws_in_core, bins=bins)
    w_hist[levels] = counts
    w_cfad[levels] = counts/(sum(counts)*0.5*1)
    counts_z, bins_z = np.histogram(zs_in_core, bins=bins_z)  
    z_hist[levels] = counts_z/(sum(counts_z)*0.5*1)
    
    
    if(len(ws_in_core) > 0):
        ninety_z[levels] = np.percentile(zs_in_core, 90)
        ninety_five_z[levels] = np.percentile(zs_in_core, 95)
        ninety_nine_z[levels] = np.percentile(zs_in_core, 99)
        ninety_w[levels] = np.percentile(ws_in_core, 90)
        ninety_five_w[levels] = np.percentile(ws_in_core, 95)
        ninety_nine_w[levels] = np.percentile(ws_in_core, 99)
    else:
        ninety_w[levels] = float('nan')
        ninety_five_w[levels] = float('nan')
        ninety_nine_w[levels] = float('nan')
        ninety_five_z[levels] = float('nan')
        ninety_nine_z[levels] = float('nan')
        ninety_z[levels] = float('nan')  

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

z_file = out_netcdf.createVariable('z', ninety_five_w.dtype, ('levels',))
z_file.long_name = 'z'
z_file.units = 'm'
z_file[:] = z_levels
                                       
out_netcdf.close()
