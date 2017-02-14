# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import math
from netCDF4 import Dataset
from scipy.interpolate import griddata
from mayavi import mlab
from datetime import datetime, timedelta
from copy import deepcopy
from matplotlib import animation
import glob

data_path = "/home/rjackson/multidop_grids/"

start_year = 2006
start_day = 19
start_month = 1
start_hour = 23
start_minute = 00
start_second = 0

end_year = 2006
end_month = 1
end_day = 20
end_hour = 1
end_minute = 0
end_second = 0

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
        print(date_str)      
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
    
# Get u,v,w,reflectivity given a time period in the CPOL dataset
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
    
    
    
    # Open multidop output
    nc = Dataset(file_name_str)
    return nc

# Get beam crossing angle - taken from DDA

def get_bca(x,y):
    berr_origin = [-12960.1,-23091.1]
    a = np.sqrt(np.multiply(x,x)+np.multiply(y,y))
    b = np.sqrt(pow(x-berr_origin[0],2)+pow(y-berr_origin[1],2))
    c = np.sqrt(berr_origin[0]*berr_origin[0]+berr_origin[1]*berr_origin[1])
    theta_1 = np.arccos(x/a)
    theta_2 = np.arccos((x-berr_origin[1])/b)
    return np.arccos((a*a+b*b-c*c)/(2*a*b))
   
def anims(time):   
    nc = get_grid_from_dda(time)
    ref = nc.variables['reflectivity'][:]
    u = nc.variables['eastward_wind'][:]
    v = nc.variables['northward_wind'][:]
    w = nc.variables['upward_air_velocity'][:]

    zh = nc.variables['z'][:]
    xh = nc.variables['x'][:]
    yh = nc.variables['y'][:]
    u = np.ma.squeeze(u)
    v = np.ma.squeeze(v)
    w = np.ma.squeeze(w)
    ref = np.squeeze(ref)
 
    grid_y, grid_z, grid_x = np.meshgrid(yh, zh, xh)

    ref = np.array(ref)
    bca = get_bca(grid_x,grid_y)
    ## Mask region outside dual doppler lobes and 
    u = np.ma.masked_where(np.logical_or(bca < math.pi/6, bca > 5*math.pi/6),u)
    no_mask = np.where(u.mask == False)
    no_mask_r = np.where(ref > 1)
    

    ## Apply mask
    grid_xm = grid_x[no_mask[0][0::10],no_mask[1][0::10],no_mask[2][0::10]]/1e3    
    grid_ym = grid_y[no_mask[0][0::10],no_mask[1][0::10],no_mask[2][0::10]]/1e3    
    grid_zm = grid_z[no_mask[0][0::10],no_mask[1][0::10],no_mask[2][0::10]]/1e3   
    um = u[no_mask[0][0::10],no_mask[1][0::10],no_mask[2][0::10]]
    vm = v[no_mask[0][0::10],no_mask[1][0::10],no_mask[2][0::10]]
    wm = w[no_mask[0][0::10],no_mask[1][0::10],no_mask[2][0::10]]
    grid_xmr = grid_x[no_mask_r[0][0::2],no_mask_r[1][0::2],no_mask_r[2][0::2]]/1e3    
    grid_ymr = grid_y[no_mask_r[0][0::2],no_mask_r[1][0::2],no_mask_r[2][0::2]]/1e3    
    grid_zmr = grid_z[no_mask_r[0][0::2],no_mask_r[1][0::2],no_mask_r[2][0::2]]/1e3   
    refm = ref[no_mask_r[0][0::2],no_mask_r[1][0::2],no_mask_r[2][0::2]]

    ## Regrid data for contour plot
    refm_gridded = griddata((grid_xmr,grid_ymr,grid_zmr), 
                            refm, 
                            (grid_x/1e3,grid_y/1e3,grid_z/1e3))
    um_gridded = griddata((grid_xm,grid_ym,grid_zm), 
                           um, 
                           (grid_x/1e3,grid_y/1e3,grid_z/1e3))
    vm_gridded = griddata((grid_xm,grid_ym,grid_zm), 
                           vm, 
                           (grid_x/1e3,grid_y/1e3,grid_z/1e3))
    wm_gridded = griddata((grid_xm,grid_ym,grid_zm), 
                           wm, 
                           (grid_x/1e3,grid_y/1e3,grid_z/1e3))

    #mlab.options.offscreen = True
    f = mlab.figure(bgcolor=(0,0,0), size=(1500, 1500))
    q = mlab.quiver3d(grid_xm, 
                      grid_ym, 
                      grid_zm, 
                      um, 
                      vm, 
                      wm,
                      extent=[-40, 40, -30, 30, 0, 20])
    
    colorbar = mlab.colorbar(q, title='m/s', orientation='vertical')
    colorbar.label_text_property.font_size = 5
    #mlab.flow(np.transpose(grid_x)/1e3, 
    #          np.transpose(grid_y)/1e3, 
    #          np.transpose(grid_z)/1e3, 
    #          np.transpose(um_gridded), 
    #          np.transpose(vm_gridded), 
    #          np.transpose(wm_gridded))

    mlab.contour3d(np.transpose(grid_x)/1e3, 
                   np.transpose(grid_y)/1e3, 
                   np.transpose(grid_z)/1e3, 
                   np.transpose(wm_gridded), colormap='blue-red',
                   contours=[-4, 4], opacity=0.5,
                   extent=[-40, 40, -30, 30, 0, 20])
    mlab.view(azimuth=72.608457142284237,
              elevation=73.373132558239334, 
              distance=201.83739701316438, 
              focalpoint=(-0.05410797, 
                          -0.34705752, 
                          10.27633847))
    mlab.title(str(time.month) + '-' + 
               str(time.day) + ' ' +
               str(time.hour) + ':' + 
               str(time.minute))
    axes = mlab.axes(extent=[-40, 40, -30, 30, 0, 20])
    axes.label_text_property.font_family = 'courier'
    axes.label_text_property.font_size = 5
    xlabel = mlab.xlabel('E of CPOL [km]')
    ylabel = mlab.ylabel('N of CPOL [km]')
    zlabel = mlab.zlabel('z [km]')
    xlabel.label_text_property.font_size = 5
    ylabel.label_text_property.font_size = 5
    zlabel.label_text_property.font_size = 5
    
    mlab.show()
    #mlab.savefig('plot3d_fig' +
    #             str(time.month) +
    #             str(time.day) +
    #             str(time.hour) + 
    #             str(time.minute) +'.png',
    #             size=(3000,2000))
    #mlab.close(f)
    


times = get_dda_times(start_year, start_month, start_day,
                      start_hour, start_minute, end_year,
                      end_month, end_day, end_hour, 
                      end_minute, minute_interval=0)

for time in times:
    anims(time)
        
#anim(times[1])
