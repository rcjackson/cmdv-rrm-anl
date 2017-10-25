import pyart
import time
import sys
import glob
import os
import time
import datetime

from dask import bag as db
from distributed import Client

radar_file_path = ('/lcrc/group/earthscience/rjackson/valentin_cpol_1b/' + 
                   sys.argv[2] + '/')
out_file_path = '/lcrc/group/earthscience/rjackson/cpol_grids_100km/'


def grid_radar(radar, grid_shape=(51, 201, 201), xlim=(-100000, 100000),
               ylim=(-100000, 100000), zlim=(0, 25000), bsp=1.0, 
               min_radius=500, h_factor=1.0, nb=1.0,
               fields=['DT', 'VT'], origin=None, gatefilter=False):
    bt = time.time()
    radar_list = [radar]
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
    grid = pyart.map.grid_from_radars(
        radar_list, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,
        weighting_function='Barnes',
        gridding_algo='map_gates_to_grid',
        h_factor=h_factor,
        min_radius=min_radius,
        bsp=bsp,
        nb=nb, 
        gatefilters=[gatefilter])
    print(time.time() - bt, 'seconds to grid radar')
    return grid


def grid_file(file_name):
    try:
        radar = pyart.io.read(file_name)
        radar_time = datetime.datetime.strptime(radar.time['units'],
            'seconds since %Y-%m-%dT%H:%M:%SZ')
        the_path = (out_file_path + '/'+ str(radar_time.year).zfill(4) +
                    '/' + str(radar_time.year).zfill(4) +
                    str(radar_time.month).zfill(2) +
                    str(radar_time.day).zfill(2) + '/')
        the_file_name = ('CPOL_GRID.' + str(radar_time.year).zfill(4) +
                         str(radar_time.month).zfill(2) +
                         str(radar_time.day).zfill(2) + '.' +
                         str(radar_time.hour).zfill(2) +
                         str(radar_time.minute).zfill(2) +
                         str(radar_time.second).zfill(2) +  '.100km.nc')
        if(os.path.isfile(the_file_name)):
            del radar
            return
        texture = pyart.retrieve.calculate_velocity_texture(radar,
            wind_size=4, vel_field='velocity')
        radar.add_field('velocity_texture', texture, replace_existing=True)
        grid = grid_radar(radar, fields=radar.fields.keys())
        del radar
        if(not os.path.isdir(the_path)):
            os.makedirs(the_path)
        pyart.io.write_grid(the_path + the_file_name, grid)
        del grid 
    except:
        print('Skipping file ' + file_name + '...corrupt.')
         
def main():
    #cluster = LocalCluster(n_workers=36, processes=True)
    #the_client = Client(cluster)
    the_client = Client(scheduler_file=sys.argv[1])
    file_list = glob.glob(radar_file_path + '/**/*.nc', recursive=True)
    the_bag = db.from_sequence(file_list)
    the_bag.map(grid_file).compute()
    the_client.shutdown()

if __name__ == '__main__':
    main()
     
