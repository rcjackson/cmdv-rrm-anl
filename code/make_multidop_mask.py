from multidop_time import do_multidop_for_time
from ipyparallel import Client
import time_procedures
from time import sleep
import time
from datetime import datetime
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import sys
import os
import math
from copy import deepcopy
from scipy import ndimage
# Get the radars from given time.
# Input the range of dates and time wanted for the collection of images

try:
    start_year = int(sys.argv[1])
    start_month = int(sys.argv[2])
    start_day = int(sys.argv[3])
    end_year = int(sys.argv[4])
    end_month = int(sys.argv[5])
    end_day = int(sys.argv[6])   
except:
   start_year = 2010
   start_day = 11
   start_month = 11
   end_year = 2011
   end_month = 5
   end_day = 15

end_hour = 17
start_hour = 0
start_minute = 1
start_second = 1
end_minute = 2
end_second = 0

# Get beam crossing angle
def get_bca(grid):
    berr_origin = [-12960.1,-23091.1]
    x,y = np.meshgrid(grid.x['data'], grid.y['data'])
    a = np.sqrt(np.multiply(x,x)+np.multiply(y,y))
    b = np.sqrt(pow(x-berr_origin[0],2)+pow(y-berr_origin[1],2))
    c = np.sqrt(berr_origin[0]*berr_origin[0]+berr_origin[1]*berr_origin[1])
    theta_1 = np.arccos(x/a)
    theta_2 = np.arccos((x-berr_origin[1])/b)
    return np.arccos((a*a+b*b-c*c)/(2*a*b))

def add_mask(grid_time):
    grid = time_procedures.get_grid_from_dda(grid_time)
    w = grid.fields['upward_air_velocity']['data']
    bca = get_bca(grid)
    le_mask = 2*np.ones(w.shape)
    # Create mask with different flags
    # 0 = include, 1 = outside of DD lobes, 2 = not in DDC
    outside_domain = np.logical_or(bca < math.pi/6, bca > 5*math.pi/6)
    num_levels = 40
    for levels in range(0,num_levels-1):
        w_outside_updraft = np.logical_or(w[levels] < 1, w[levels] > 99.0)
        le_mask[levels,outside_domain] = 1
        w[levels] = np.ma.masked_where(np.logical_or(w_outside_updraft,
                                                     outside_domain), w[levels])
        
    grid_z = grid.point_z['data']
    # Set mask to exclude data outside of updrafts
    in_updraft = np.where(np.logical_and(bca > math.pi/6, bca < 5*math.pi/6))
    w_temp = deepcopy(w)
    w_temp[~w.mask] = 1
    w_temp[w.mask] = 0
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
    
    # Find deep convective cores 
    for labels in range(1, len(max_z)-1):
        indicies = np.ma.where(updrafts == labels)                                
                        
        if(len(indicies[0]) > 0 and 
           max_z[labels] >= 15000 and 
           min_z[labels] <= 1000):
            le_mask[updrafts == labels] = 0
    
    field_mask = {'data': le_mask, 
                  'units': '0 = include, 1 = not in lobes, 2 = not in deep core'}        
    grid.add_field('mask', field_mask, replace_existing=True)
    time_procedures.write_grid(grid_time, grid)

    
        
times = time_procedures.get_dda_times(start_year, start_month, start_day,
                                      start_hour, start_minute, end_year,
                                      end_month, end_day, end_minute, 
                                      end_second)
  
print(times)
print('About to process ' + str(len(times)) + 'grids')
serial = 1 

if(serial == 0):
    # Get iPython cluster
    state = 0
    while state == 0:
        try:
            My_Cluster = Client()
            My_View = My_Cluster[:]
            state = 1
        except:
            state = 0
            print('Cluster not ready for me')
            sleep(10)

    #Turn off blocking so all engines can work async
    My_View.block = False

    #on all engines do an import of Py-ART
    My_View.execute('import matplotlib')
    My_View.execute('matplotlib.use("agg")')
    My_View.execute('import pyart')
    My_View.execute('import numpy as np')
    t1 = time.time()

    #for timer in times:
    #  do_multidop_for_time(timer)

    #Map the code and input to all workers
    result = My_View.map_async(add_mask, times)

    #Reduce the result to get a list of output
    qvps = result.get()
    tt = time.time() - t1
    print(tt)
    print(tt/len(times))
else:
    for timer in times:
        add_mask(timer)


    
