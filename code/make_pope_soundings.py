#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 15:28:19 2016

@author: rjackson
"""

import sklearn
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta

# Load sounding file
sounding_file = '/home/rjackson/data/soundings/twpsondewnpnC3.b1.20020401.043300..20091230.172000.custom.cdf'

Sounding_netcdf = Dataset(sounding_file, mode='r')

# Convert timestamps to datetime format
Time = Sounding_netcdf.variables['time_offset'][:]
base_time = Sounding_netcdf.variables['base_time'][:]
p = Sounding_netcdf.variables['pres'][:]
u = Sounding_netcdf.variables['u_wind'][:]
v = Sounding_netcdf.variables['v_wind'][:]
t = Sounding_netcdf.variables['tdry'][:]

base_timestamp = timedelta(seconds=float(base_time)) + datetime(1970,1,1)

# Restrict analysis to September to April and 2350 UTC.

start_year = 2002
end_year = 2009
pres_levels = [1013, 950, 925, 900, 850, 800, 750, 700, 
               650, 600, 550, 400, 300, 200, 100]

def find_nearest(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

## Save soundings at 16 levels for later

i = 0
time_stamps = []
u_soundings = []
v_soundings = []
t_soundings = []
pres_index = []
time_index = 0

while(time_index < len(Time)):
    time_stamp = timedelta(seconds=float(Time[time_index])) + base_timestamp

    if(time_stamp.year >= start_year and time_stamp.year <= end_year):
        if(time_stamp.month <= 4 or time_stamp.month >= 9):
            if(time_stamp.hour == 23):
                
                nearest_index = time_index
                furthest_index = time_index + 1
                
                while(p[nearest_index-1] > p[nearest_index]):
                    nearest_index = nearest_index - 1  
                # Go forwards in time until altitude_stops_increasing
                
                while(p[furthest_index+1] < p[furthest_index]):
                    furthest_index = furthest_index + 1  
                time_stamps.append(datetime(time_stamp.year, time_stamp.month, 
                                           time_stamp.day, 23, 50, 00))
                
                us = u[int(nearest_index):int(furthest_index)]
                vs = v[int(nearest_index):int(furthest_index)]
                pressure = p[int(nearest_index):int(furthest_index)]
                ts = t[int(nearest_index):int(furthest_index)]
            
                # Limit to 16 pressure levels
            
                for j in range(0, len(pres_levels)-1):
                    pres_index.append([pres_index, find_nearest(pressure, pres_levels[j])])
                
                
                year_str = "%04d" % time_stamp.year
                month_str = "%02d" % time_stamp.month
                day_str = "%02d" % time_stamp.day
                sounding_file_name = ('/home/rjackson/data/soundings/pope_soundings/' +
                                      year_str + month_str + day_str)
                print(sounding_file_name)
                file = open(sounding_file_name, 'w')

                # Take levels from the sounding and place them into the file
                us = us[~us.mask]
                vs = vs[~vs.mask]
                pressure = pressure[~us.mask]
 
                for i in pres_index:
                    input_string = (str(pressure[i]) + ' ' + str(u[i]) + ' ' + str(v[i]) + '\n')
                    file.write(input_string)

                file.close()
                time_index = furthest_index + 1
            else:
                time_index = time_index + 1
        else:
            time_index = time_index + 1
    else:
        time_index = time_index + 1