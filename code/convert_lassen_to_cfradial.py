## Converts lassen files in given directory to Cf/Radial format

import matplotlib
# Needed for blues
matplotlib.use('agg')
import pyart
import glob
import os
import time
from datetime import datetime
from parse import *
from ipyparallel import Client
from time import sleep

data_file_path = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_lassen/'
out_file_path = '/home/rjackson/data/radar/cpol/'

# Get list of radar files (recursively)
file_list = glob.glob(data_file_path + '/**/*_PPI.lassen', recursive=True)
print('Converting ' + str(len(file_list)) + ' files to Cf/Radial')

# Convert each radar file
def convert_file(radar_file):
    import os
    from parse import parse
    out_file_path = '/home/rjackson/data/radar/cpol'
    print('Reading ' + radar_file)
    Radar = pyart.io.read_rsl(radar_file)
    Radar.info()
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

    if(not os.path.exists(out_file_path_this)):
        os.makedirs(out_file_path_this)
    
    out_file_name = (out_file_path_this +
                     'Gunn_pt_' +
                     year_str +
                     month_str +
                     day_str +
                     hour_str +
                     minute_str +
                     second_str +
                     '.nc')
    print('Writing ' + out_file_name)
    pyart.io.write_cfradial(out_file_name, Radar)


serial = 0

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

    #Map the code and input to all workers
    result = My_View.map_async(convert_file, file_list)

    #Reduce the result to get a list of output
    qvps = result.get()
    tt = time.time() - t1
    print(tt)
    print(tt/len(times))
else:
    for file_name in file_list:
        convert_file(file_name)








