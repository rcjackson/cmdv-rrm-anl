#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 11:37:23 2017

@author: rjackson
"""
from matplotlib import use
use('agg')
import pyart
from netCDF4 import Dataset
import numpy as np
import math
import os.path
from glob import glob
from datetime import datetime, timedelta
from dask import bag as db
from distributed import Client, LocalCluster

map_path = '100km_cfradial_to_scrib_mapping.nc'
echo_top_data_path = '/lcrc/group/earthscience/rjackson/echo_tops/'
century_weight_file = '100km_weights.nc'
out_path = '/lcrc/group/earthscience/rjackson/echo_top_rrm_grid/'

def convert_echo_tops(echo_top_data_file):
    echo_top_dataset = Dataset(echo_top_data_file)
    basetime = echo_top_dataset['time'].units
    basetime = datetime.strptime(basetime, 'seconds since %Y-%m-%d %H:%M:%S')
    the_shape = echo_top_dataset['cpol_T'][:].shape
    for i in range(the_shape[0]):
        ctop = echo_top_dataset['cpol_T'][i]
        time_after = timedelta(seconds=float(echo_top_dataset['time'][i])) + basetime
        out_file_name = (out_path + 'cloudtop_twprrm' + 
                         time_after.strftime('%Y%m%d.%H%M%S') + '.nc')
        if(not os.path.isfile(out_file_name)):
            print('Processing scan ' + str(i))
            ctop_dest = np.nan*np.zeros(dest_centers_lat.shape)
            count =  np.zeros(dest_centers_lat.shape)
            equivalent_point = [np.where(points == col[i]-1) for i in range(0, len(S))]
            for i in range(0, len(S)):
                if(len(equivalent_point[i][0]) == 1):
                    ctop_src = ctop[equivalent_point[i][0], equivalent_point[i][1]]
                    if(ctop_src > -999.0):
                        if(not np.isfinite(ctop_dest[row[i]-1])):
                            ctop_dest[row[i]-1] = 0

                    ctop_dest[row[i]-1] = ctop_dest[row[i]-1] + S[i]*ctop_src
                    count[row[i]-1] += 1

            ctop_dest[frac_b != 0] = ctop_dest[frac_b != 0]/frac_b[frac_b != 0] 
            ctop_dest[count < 50] = np.nan 

            # Create an output SCRIP file
            new_dataset = Dataset(out_file_name, mode='w')
            new_dataset.createDimension('grid_size', grid_size)
            new_dataset.createDimension('grid_corners', num_corners)
            new_dataset.createDimension('grid_rank', grid_rank)

            radar_est_cloud_top = new_dataset.createVariable('radar_est_cloud_top',
                                                             ctop_dest.dtype,
                                                            'grid_size')
            radar_est_cloud_top.units = 'm'
            radar_est_cloud_top.long_name = 'Radar estimated cloud top height'
            radar_est_cloud_top[:] = ctop_dest

            grid_center_lat = new_dataset.createVariable('grid_center_lat',
                                           dest_corners_lat.dtype,
	                                  'grid_size')
            grid_center_lat.units = 'degrees'
            grid_center_lat[:] = dest_centers_lat

            grid_center_lon = new_dataset.createVariable('grid_center_lon',
                                               dest_corners_lat.dtype,
                                              'grid_size')
            grid_center_lat.units = 'degrees'
            grid_center_lon[:] = dest_centers_lon

            grid_corner_lat = new_dataset.createVariable('grid_corner_lat',
                                            dest_corners_lat.dtype,
                                            ('grid_size',
                                            'grid_corners'),
                                            fill_value=float(9.97e36))
            grid_corner_lat.units = 'degrees'
            grid_corner_lat[:] = dest_corners_lat

            grid_corner_lon = new_dataset.createVariable('grid_corner_lon',
                                                         dest_corners_lat.dtype,
                                                        ('grid_size',
                                                        'grid_corners'),
                                                        fill_value=float(9.97e36))
            grid_corner_lon.units = 'degrees'
            grid_corner_lon[:] = dest_corners_lon
            grid_dims = new_dataset.createVariable('grid_dims', float, 'grid_rank')
            grid_dims[:] = 1
            grid_imask = new_dataset.createVariable('grid_imask',
	                                            dest_corners_lat.dtype,
                                                    'grid_size', 
                                                    fill_value=float(9.97e36))
            grid_imask[:] = dest_mask
            new_dataset.close()
        else:
            print(out_file_name + ' Already exists...skipping!')

         
if __name__ == '__main__':
    century_weights_dataset = Dataset(century_weight_file)
    map_dataset = Dataset(map_path)
        
    points = map_dataset['scrib_index'][:]
    dest_corners_lat = century_weights_dataset['yv_b'][:]
    dest_centers_lat = century_weights_dataset['yc_b'][:]
    dest_corners_lon = century_weights_dataset['xv_b'][:]
    dest_centers_lon = century_weights_dataset['xc_b'][:]
    dest_mask = century_weights_dataset['mask_b'][:]
    dest_area = century_weights_dataset['area_b'][:]
    S = century_weights_dataset['S'][:]
    row = century_weights_dataset['row'][:]
    col = century_weights_dataset['col'][:]
    frac_b = century_weights_dataset['frac_b'][:]
    num_corners = dest_corners_lat.shape[1]
    grid_size = dest_corners_lat.shape[0]
    grid_rank = 1
    grid_dims = 1
    file_list = glob(echo_top_data_path + '*.cdf')
    the_bag = db.from_sequence(file_list)
    the_bag.map(convert_echo_tops).compute()
        
