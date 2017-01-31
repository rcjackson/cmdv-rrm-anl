## Converts lassen files in given directory to Cf/Radial format

import matplotlib
# Needed for blues
matplotlib.use('agg')
import pyart
import glob
from datetime import datetime
from parse import *

data_file_path = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_lassen/'
out_file_path = '/home/rjackson/data/radar/cpol/'

# Get list of radar files (recursively)
file_list = glob.glob(data_file_path + '/**/*_PPI.lassen', recursive=True)
print('Converting ' + str(len(file_list)) + ' files to Cf/Radial')

# Convert each radar file
for radar_file in file_list:
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
    out_file_name = (out_file_path + 
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


  






