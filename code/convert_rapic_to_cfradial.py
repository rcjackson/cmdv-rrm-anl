## Converts rapic files in given directory to Cf/Radial format

import pyart
import glob

data_file_path = '/home/rjackson/data/radar/berr/berr_2009'

# Get list of radar files (recursively)

file_list = glob.glob(data_file_path + '/**/*.rapic', recursive=True)
print('Converting ' + str(len(file_list)) + ' files to Cf/Radial')

for radar_file in file_list:
    Radar = pyart.io.read_rsl(radar_file)
    out_file_name = radar_file[0:-5] + '.cfradial'
    pyart.io.write_cfradial(out_file_name, Radar)


  






