import pyart
import numpy as np
import glob
import dask.bag as db
import xarray
import time
import sys
import pandas

from distributed import Client, LocalCluster
from dask import delayed
from datetime import datetime
from netCDF4 import Dataset


def get_rain_rate_over_spots(file, disdrometer_index, 
                             profiler_index):
    
    grid = pyart.io.read_grid(file)
    rain_dis = grid.fields['radar_estimated_rain_rate']['data'][3,
        disdrometer_index[0], disdrometer_index[1]]
    rain_prof = grid.fields['radar_estimated_rain_rate']['data'][3, 
        profiler_index[0], profiler_index[1]]
    echo = grid.fields['radar_echo_classification']['data'][:,
        profiler_index[0], profiler_index[1]]
    if(np.any(echo == 7) or np.any(echo == 8)):
        graupel = 1
    else:
        graupel = 0
    if(np.any(echo == 9)):
        hail = 1
    else:
        hail = 0
    rain_dis = np.squeeze(rain_dis)
    rain_prof = np.squeeze(rain_prof)
    print(file + ' successful!')
    del grid
    #except:
    #    print(file + 'is corrupt...skipping!')
    #    return np.nan*np.ones(4)
                       

    return np.array([rain_dis, rain_prof, graupel, hail])

def parse_time(file):
    return datetime.strptime(file[-28:-15], '%Y%m%d_%H%M')

def ten_minute_average_impact(inp_time):
    disdrometer_path = '/lcrc/group/earthscience/rjackson/disdrometer/'
    file_name = ('twpdisdrometerC3.b1.' + "%04d" % inp_time.year +
                 "%02d" % inp_time.month + "%02d" % inp_time.day +
                 ".*.cdf")
    try:
        ds = xarray.open_mfdataset(disdrometer_path+file_name)
        ds.resample(time='10min')
        rain_value = ds['rain_rate'].sel(time=inp_time).values
        print(str(inp_time) + ' sucessfully loaded!')
        del ds
    except (IOError, KeyError, ValueError):
        return np.nan
    return rain_value  


def ten_minute_average_vdis(inp_time):
    disdrometer_path = '/lcrc/group/earthscience/rjackson/disdrometer/'
    file_name = ('twpvdisC3.b1.' + "%04d" % inp_time.year +
                 "%02d" % inp_time.month + "%02d" % inp_time.day +
                 ".*.cdf")
    try:
        ds = xarray.open_mfdataset(disdrometer_path+file_name)
        ds.resample(time='10min')
        rain_value = ds['rain_rate'].sel(time=inp_time).values
        print(str(inp_time) + ' sucessfully loaded!')
        del ds
    except (IOError, KeyError, ValueError):
        return np.nan
    return rain_value

    
def get_dros_class(inp_time):
    file_path = '/home/rjackson/data/Drosdowsky.cdf'
    in_netcdf = Dataset(file_path)
    year = in_netcdf.variables['year'][:]
    month = in_netcdf.variables['month'][:]
    day = in_netcdf.variables['day'][:]
    groups = in_netcdf.variables['groups'][:]
    the_index = np.where(np.logical_and.reduce((
        year == inp_time.year, month == inp_time.month, day == inp_time.day)))
    if(the_index[0]):
        return groups[the_index[0][0]]
    else:
        return np.nan

    
def get_mjo_index(inp_time):
    mjo_index_file = '/home/rjackson/data/rmm.74toRealtime.txt'
    data = pandas.read_csv(mjo_index_file,
                           header=2,
                           delim_whitespace=True)
    data_matrix = np.ma.array(data.values)
    yearm = data_matrix[:, 0]
    monthm = data_matrix[:, 1]
    daym = data_matrix[:, 2]
    index = data_matrix[:, 5]
    the_index = np.where(np.logical_and.reduce((
        yearm == inp_time.year, monthm == inp_time.month, daym == inp_time.day)))
    if(the_index[0]):
        return index[the_index[0][0]]
    else:
        return np.nan

# Get information about disdrometer and profiler indicies in advance to speed
# up algorithm
if __name__ =='__main__':
    year = sys.argv[1]
    file_path = ('/lcrc/group/earthscience/radar/CPOL_level_1b/' +
                 'GRIDDED/GRID_70km_1000m/' + year + '/')
    scheduler_file = '/home/rjackson/scheduler.json'
    dis_lat = -12.0 - 25/60.0 - 30/3600.0
    dis_lon = 130.0 + 53.0/60.0 + 31.2/3600.0
    profiler_lat = dis_lat
    profiler_lon = dis_lon
    #client = Client(scheduler_file=scheduler_file)
    #cluster = LocalCluster(n_workers=36)
    #client = Client(cluster)
    file_list = glob.glob(file_path + '/**/*.nc', recursive=True)
    first_grid = pyart.io.read_grid(file_list[0])
    lats = first_grid.point_latitude['data'][0,:,0]
    lons = first_grid.point_longitude['data'][0,0,:]
    lat_index_dis = 50
    lon_index_dis = 53
    print(lat_index_dis, lon_index_dis)
    profiler_lat_ind= (np.abs(lats-profiler_lat)).argmin()
    profiler_lon_ind = (np.abs(lons-profiler_lon)).argmin()
    dis_ind = (lat_index_dis,lon_index_dis)
    prof_ind = (profiler_lat_ind, profiler_lon_ind)
    
    bt = time.time()
    print('## Parsing dates and times...')
    file_list = sorted(file_list)
    file_bag = db.from_sequence(file_list)
    time_list = file_bag.map(parse_time).compute()
    print('## Getting Drosdowsky classification...')
    time_bag = db.from_sequence(time_list)
    dros_class = time_bag.map(get_dros_class).compute()
    dros_class = np.array(dros_class)
    
    mjo_index = time_bag.map(get_mjo_index).compute()
    mjo_index = np.array(mjo_index)
    print('## First, aggregating the video disdrometer observations to 10 min')   
    vdis_rainfall = time_bag.map(ten_minute_average_vdis).compute()
    vdis_rainfall = np.stack(vdis_rainfall)
    print('## Aggregating impact disdrometer data')
    dis_rainfall = time_bag.map(ten_minute_average_impact).compute()
    dis_rainfall = np.stack(dis_rainfall)
    
    print('## Getting rainfall values...')
    dis_values = file_bag.map(lambda x: get_rain_rate_over_spots(
                              x, dis_ind, dis_ind))
    dis_values = np.stack(dis_values, axis=0)
    print(dis_values.shape)
    dis_rain = dis_values[:,0]
    prof_rain = dis_values[:,1]
    graupel = dis_values[:,2]
    hail = dis_values[:,3]

    ds = xarray.Dataset({'rainfall_cpol_over_impact': (['time'], dis_rain),
                         'rainfall_cpol_over_vdis': (['time'], prof_rain),
                         'rainfall_vdis': ('time', vdis_rainfall),
                         'rainfall_disdro': ('time', dis_rainfall),
                         'graupel_present': ('time', graupel),
                         'hail_present': ('time', hail),
                         'dros_class': ('time', dros_class),
                         'mjo_index': ('time', mjo_index)},
                         coords={'time': time_list},
                         attrs={'units': 'mm hr-1', 'long_name':
                                ('Rainfall rate algorithm based on Thompson' +
                                 'et al. 2016.')})
    print(ds)
    ds.to_netcdf(path=('rainfall_rate_timeseries' + year + '.cdf'), mode='w')
    print('## Task completed in ' + str((time.time-bt)/60.0))
    client.shutdown()
