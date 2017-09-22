import pyart
import numpy as np
import glob
import dask.bag as db
import xarray
import time

from distributed import Client, LocalCluster
from dask import delayed
from datetime import datetime


def get_rain_rate_over_spots(file, disdrometer_index, 
                             profiler_index):
    try:
        grid = pyart.io.read_grid(file)
        rain_dis = grid.fields['radar_estimated_rain_rate']['data'][:,
            disdrometer_index[0], disdrometer_index[1]]
        rain_prof = grid.fields['radar_estimated_rain_rate']['data'][:, 
            profiler_index[0], profiler_index[1]]
        rain_dis = np.squeeze(rain_dis)
        rain_prof = np.squeeze(rain_prof)
        del grid
    except:
        print(file + 'is corrupt...skipping!')
        return np.nan*np.ones((2,41))
                       

    return np.array([rain_dis, rain_prof])

def parse_time(file):
    return datetime.strptime(file[-28:-15], '%Y%m%d_%H%M')


# Get information about disdrometer and profiler indicies in advance to speed
# up algorithm
if __name__ =='__main__':
    file_path = ('/lcrc/group/earthscience/rjackson/valentin_cpol_1b/' +
                 'GRIDDED/GRID_70km_1000m/')
    scheduler_file = '/home/rjackson/scheduler.json'
    dis_lat = -12.0 - 25/60.0 - 30/3600.0
    dis_lon = 130.0 + 53.0/60.0 + 31.2/3600.0
    profiler_lat = -12.444857
    profiler_lon = 130.957718
    client = Client(scheduler_file=scheduler_file)
    #cluster = LocalCluster(n_workers=36)
    #client = Client(cluster)

    file_list = glob.glob(file_path + '/**/*.nc', recursive=True)
    first_grid = pyart.io.read_grid(file_list[0])
    lats = first_grid.point_latitude['data'][0,:,0]
    lons = first_grid.point_longitude['data'][0,0,:]
    lat_index_dis = (np.abs(lats-dis_lat)).argmin()
    lon_index_dis = (np.abs(lons-dis_lon)).argmin()
    profiler_lat_ind= (np.abs(lats-profiler_lat)).argmin()
    profiler_lon_ind = (np.abs(lons-profiler_lon)).argmin()
    dis_ind = (lat_index_dis,lon_index_dis)
    prof_ind = (profiler_lat_ind, profiler_lon_ind)
    
    bt = time.time()
    print('## Parsing dates and times...')
    file_bag = db.from_sequence(file_list)
    time_list = file_bag.map(parse_time).compute()


    print('## Getting rainfall values...')
    dis_values = file_bag.map(lambda x: get_rain_rate_over_spots(
                              x, dis_ind, prof_ind))
    dis_values = np.stack(dis_values, axis=0)
    dis_rain = dis_values[:,0,:]
    prof_rain = dis_values[:,1,:]

    ds = xarray.Dataset({'rainfall_disdrometer': (['time', 'z'], dis_rain),
                         'rainfall_profiler': (['time', 'z'], prof_rain)},
                         coords={'height': (['z'], first_grid.z['data']),
                                 'time': time_list},
                         attrs={'units': 'mm hr-1', 'long_name':
                                ('Rainfall rate algorithm based on Thompson' +
                                 'et al. 2016.')})
    print(ds)
    ds.to_netcdf(path='rainfall_rates_disdrometer.cdf', mode='w')
    print('## Task completed in ' + str((time.time-bt)/60.0))
    client.shutdown()
