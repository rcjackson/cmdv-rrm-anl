from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
from copy import deepcopy
import math
import dask.array as da
from distributed import Client, LocalCluster
from dask import delayed, compute
import time
import sys
from scipy import ndimage
import pandas
import time_procedures
import matplotlib
matplotlib.use('Agg')
import pyart

# Start a cluster with x workers
cluster = LocalCluster(n_workers=int(sys.argv[1]), processes=False)
client = Client(cluster)

# Input the range of dates and time wanted for the collection of images
start_year = 2005
start_month = 11
start_day = 1
start_hour = 1
start_minute = 0
start_second = 0

end_year = 2011
end_month = 5
end_day = 2
end_hour = 0
end_minute = 00
end_second = 0

data_path = '/lcrc/group/earthscience/rjackson/multidop_grids/ddop/'


# Get beam crossing angle between radars
def get_bca(grid):
    berr_origin = [-12960.1, -23091.1]
    x, y = np.meshgrid(grid.x['data'], grid.y['data'])
    a = np.sqrt(np.square(x) + np.square(y))
    b = np.sqrt(pow(x - berr_origin[0], 2)+pow(y - berr_origin[1], 2))
    c = np.sqrt(berr_origin[0]*berr_origin[0] + berr_origin[1]*berr_origin[1])
    theta_1 = np.arccos(x/a)
    theta_2 = np.arccos((x - berr_origin[1])/b)
    return np.arccos((a*a + b*b - c*c)/(2*a*b))


def get_updrafts(time):
    pyart_grid = time_procedures.get_grid_from_dda(time)
    w = pyart_grid.fields['upward_air_velocity']['data']

    for levels in range(0, num_levels-1):
        w_outside_updraft = np.logical_or(w[levels] < 1, w[levels] > 99.0)
        outside_dd_lobes = np.logical_or(bca < math.pi/6, bca > 5*math.pi/6)
        w[levels] = np.ma.masked_where(
                    np.logical_or(w_outside_updraft,
                                  outside_dd_lobes), w[levels])

    grid_z = pyart_grid.point_z['data']

    # Set mask to exclude data outside of updrafts
    w_temp = deepcopy(w)
    w_temp[~w_temp.mask] = 1
    w_temp[w_temp.mask] = 0
    w_temp.mask = False

    six_connected_structure = [[[0, 0, 0],
                                [0, 1, 0],
                                [0, 0, 0]],
                               [[0, 1, 0],
                                [1, 1, 1],
                                [0, 1, 0]],
                               [[0, 0, 0],
                                [0, 1, 0],
                                [0, 0, 0]]]
    updrafts, num_updrafts = ndimage.measurements.label(
                             w_temp, structure=six_connected_structure)

    # Get statistics in continous regions
    index = np.arange(0, num_updrafts + 1)
    max_z = ndimage.measurements.maximum(grid_z,
                                         labels=updrafts,
                                         index=index)
    min_z = ndimage.measurements.minimum(grid_z,
                                         labels=updrafts,
                                         index=index)

    max_w_individual = []
    level_individual = []

    # Find deep convective cores and get max updraft speeds
    for levels in range(0, num_levels-1):
        label_level = updrafts[levels]
        masked_array = np.ma.zeros(updrafts.shape)
        masked_array.mask = True
        w_temp = w[levels]

        for labels in range(1, len(max_z)-1):
            indicies = np.ma.where(label_level == labels)

            if(len(indicies[0]) > 0 and
               max_z[labels] >= 15000 and
               min_z[labels] <= 1000):
                max_w_individual.append(max(w_temp[indicies]))
                level_individual.append(levels)

    # Convert to list of individual max w's for each updraft
    max_w_individual = np.array(max_w_individual)
    level_individual = np.array(level_individual)

    return_array = np.ma.zeros((len(max_w_individual), 3))
    return_array[:, 0] = max_w_individual
    return_array[:, 1] = level_individual
    return return_array


# Get the radars for a specific time
times = time_procedures.get_dda_times(start_year, start_month, start_day,
                                      start_hour, start_minute, end_year,
                                      end_month, end_day, end_hour,
                                      end_minute, minute_interval=0)

# Load Pope regimes
in_netcdf = Dataset('/home/rjackson/data/Drosdowsky.cdf',
                    mode='r')
year = in_netcdf.variables['year'][:]
month = in_netcdf.variables['month'][:]
day = in_netcdf.variables['day'][:]
groups = in_netcdf.variables['groups'][:]

drosdates = []
for i in range(0, len(day)):
    drosdates.append(datetime(year=int(year[i]),
                              month=int(month[i]),
                              day=int(day[i])))

# Since grids are uniform, calculate beam crossing angle for first grid and
# apply to all
first_grid = time_procedures.get_grid_from_dda(times[0])
bca = get_bca(first_grid)
num_levels = 40
z_levels = first_grid.z['data']
count = 0
dros_regime = int(sys.argv[2])

# Filter out data not in Pope regime
dros_times = []
for time in times:
    # Look for date in Pope regime data
    cur_date = datetime(year=time.year, month=time.month, day=time.day)
    inds = np.where([day <= cur_date for day in drosdates])
    dros_index = inds[0][-1]
    if(groups[dros_index] == dros_regime):
        print((drosdates[dros_index], time))
        dros_times.append(time)

in_netcdf.close()

# Get delayed structure to load files in parallel
get_file = delayed(get_updrafts)

# Calculate PDF
mean_w = np.ma.zeros(num_levels)
median_w = np.ma.zeros(num_levels)
ninety_w = np.ma.zeros(num_levels)
ninety_five_w = np.ma.zeros(num_levels)
ninety_nine_w = np.ma.zeros(num_levels)
print('Doing parallel grid loading...')
t1 = time.time()
ws = []

for i in range(0, len(dros_times), int(len(dros_times)/4)):
    ws_temp = [get_file(times)
               for times in dros_times[i:(i+int(len(dros_times)/4))]]
    ws_temp = compute(*ws_temp)
    ws.append(ws_temp)

# Print chunk sizes
for arrays in ws:
    array_temp = np.concatenate(arrays)
    print(array_temp.shape)

ws = np.concatenate([np.concatenate(arrays) for arrays in ws])
t2 = time.time() - t1
print('Total time in s: ' + str(t2))
print('Time per scan = ' + str(t2/len(dros_times)))
level_individual = ws[:, 1]
w_individual = ws[:, 0]

for levels in range(0, num_levels):
    w_new = w_individual[level_individual == levels]
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

print(mean_w)
print('Writing netCDF file...')

# Save to netCDF file
out_netcdf = Dataset('wpdfdros' + str(dros_regime) + '_.cdf', 'w')
out_netcdf.createDimension('levels', num_levels)
mean_file = out_netcdf.createVariable(
            'mean', mean_w.dtype, ('levels',))
mean_file.long_name = 'Mean w'
mean_file.units = 'm s-1'
mean_file[:] = mean_w

median_file = out_netcdf.createVariable(
              'median', median_w.dtype, ('levels',))
median_file.long_name = 'median w'
median_file.units = 'm s-1'
median_file[:] = median_w

ninety_file = out_netcdf.createVariable(
              'ninety', ninety_w.dtype, ('levels',))
ninety_file.long_name = '90% w'
ninety_file.units = 'm s-1'
ninety_file[:] = ninety_w

n5_file = out_netcdf.createVariable(
          'ninety_five', ninety_five_w.dtype, ('levels',))
n5_file.long_name = '95W w'
n5_file.units = 'm s-1'
n5_file[:] = ninety_five_w

n5_file = out_netcdf.createVariable(
          'ninety_nine', ninety_five_w.dtype, ('levels',))
n5_file.long_name = '99W w'
n5_file.units = 'm s-1'
n5_file[:] = ninety_nine_w

z_file = out_netcdf.createVariable(
         'z', ninety_five_w.dtype, ('levels',))
z_file.long_name = 'z'
z_file.units = 'm'
z_file[:] = z_levels

out_netcdf.close()
