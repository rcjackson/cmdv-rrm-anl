## Converts lassen files in given directory to Cf/Radial format

import matplotlib
# Needed for blues
matplotlib.use('agg')
import pyart
import glob
import fnmatch
import os
import time
from datetime import datetime
from parse import *
from ipyparallel import Client
from time import sleep

data_file_path = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_lassen/'
out_file_path = '/lcrc/group/earthscience/rjackson/cpol/'

list_file = open('list_of_files0203')

# Get list of radar files (recursively) (Python 2.7)
file_list = list_file.readlines() 
list_file.close()

# Python 3.5 only (uncomment if using)
#file_list = glob.glob(data_file_path + '/**/*_PPI.lassen', recursive=True)
print('Converting ' + str(len(file_list)) + ' files to Cf/Radial')
print(file_list[1])
# Convert each radar file
def convert_file(radar_file):
    import os
    from parse import parse
    out_file_path = '/lcrc/group/earthscience/rjackson/cpol/'
    
    print('Reading ' + radar_file)
    try:
        Radar = pyart.io.read_rsl(radar_file)
        time_string = Radar.time['units']
        parse_string = 'seconds since {:ti}'
        radar_datetime = parse(parse_string, time_string)
        radar_datetime = radar_datetime[0]

        year_str = "%04d" % radar_datetime.year 
        month_str = "%02d" % radar_datetime.month
        day_str = "%02d" % radar_datetime.day
        hour_str = "%02d" % radar_datetime.hour
        minute_str = "%02d" % radar_datetime.minute
        second_str = "%02d" % radar_datetime.second
      
        out_file_path_this = (out_file_path +
                              '/' +
                              year_str +
                              '/' +
                              year_str +
                              month_str +
                              day_str +
                              '/')

        out_file_name = (out_file_path_this +
                        'Gunn_pt_' +
                         year_str +
                         month_str +
                         day_str +
                         hour_str +
                         minute_str +
                         second_str +
                         Radar.scan_type +
                         '.nc')

        if(not os.path.isfile(out_file_name)):
            if(not os.path.exists(out_file_path_this)):
                os.makedirs(out_file_path_this)
    
            print('Writing ' + out_file_name)
            pyart.io.write_cfradial(out_file_name, Radar)
        else:
            print('Skipping ' + radar_file + ', converted output already exists.')
    except:
        print('Skipping ' + radar_file + ' cannot be read by RSL.')

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
    My_View.execute('import time_procedures')
    t1 = time.time()

    #for rad_date in dates:
    #    display_time(rad_date)
    print('Checking for already converted files...')
    file_list1 = []
    for file_name in file_list:
        print(file_name[-4:])
        if(file_name[-4:] == 'vol'):
            hex_timestamp = file_name[-13:-5]
            print(hex_timestamp)
            time_stamp = datetime.utcfromtimestamp(float(int(hex_timestamp, 16)))     
            year_str = "%04d" % time_stamp.year
            month_str = "%02d" % time_stamp.month
            day_str = "%02d" % time_stamp.day
            hour_str = "%02d" % time_stamp.hour
            minute_str = "%02d" % time_stamp.minute
            second_str = "%02d" % time_stamp.second
        else:
            date_str = file_name[-26:-12]
            year_str = date_str[0:4]
            month_str = date_str[4:6]
            day_str = date_str[6:8]
            hour_str = date_str[8:10]
            minute_str = date_str[10:12]
            second_str = date_str[12:14]
       
        out_file_path_this = (out_file_path +
                              '/' +
                              year_str +
                              '/' +
                              year_str +
                              month_str +
                              day_str +
                              '/')

        out_file_name = (out_file_path_this +
                         'Gunn_pt_' +
                          year_str +
                          month_str +
                          day_str +
                          hour_str +
                          minute_str +
                          second_str +
                         'ppi.nc')
        if(not os.path.isfile(out_file_name)):
            file_list1.append(file_name[:-1])

    print('Starting conversion...')
    #Map the code and input to all workers
    result = My_View.map_async(convert_file, file_list1)

    #Reduce the result to get a list of output
    qvps = result.get()
    tt = time.time() - t1
    print(tt)
    print(tt/len(times))
else:
    for file_name in file_list:
        if(file_name[-4:-1] == 'vol'):
            hex_timestamp = file_name[-13:-5]
            print(hex_timestamp)
            time_stamp = datetime.utcfromtimestamp(float(int(hex_timestamp, 16)))     
            year_str = "%04d" % time_stamp.year
            month_str = "%02d" % time_stamp.month
            day_str = "%02d" % time_stamp.day
            hour_str = "%02d" % time_stamp.hour
            minute_str = "%02d" % time_stamp.minute
            second_str = "%02d" % time_stamp.second
        else:
            date_str = file_name[-26:-12]
            year_str = date_str[0:4]
            month_str = date_str[4:6]
            day_str = date_str[6:8]
            hour_str = date_str[8:10]
            minute_str = date_str[10:12]
            second_str = date_str[12:14]
       
        out_file_path_this = (out_file_path +
                              '/' +
                              year_str +
                              '/' +
                              year_str +
                              month_str +
                              day_str +
                              '/')

        out_file_name = (out_file_path_this +
                         'Gunn_pt_' +
                          year_str +
                          month_str +
                          day_str +
                          hour_str +
                          minute_str +
                          second_str +
                         'ppi.nc')
        
        if(not os.path.isfile(out_file_name)):
            convert_file(file_name[:-1])
        else:
            print('Skipping ' + year_str + '-' + day_str + '-' + month_str)








