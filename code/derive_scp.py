import matplotlib
matplotlib.use('Agg')
import pyart
from matplotlib import pyplot as plt
import numpy as np
import glob
import os
from copy import deepcopy
from ipyparallel import Client
from time import sleep
import time
import time_procedures
from netCDF4 import Dataset

## Berrima - 2009-2011 (new format), 2005-2005 (old format)
start_year = 2011
start_month = 1
start_day = 28
start_hour = 0
start_minute = 1

end_year = 2011
end_month = 1
end_day = 29
end_hour = 1
end_minute = 2

def SCP(rad_time):
    import pyart
    import matplotlib
    import os
    from scipy import ndimage
    os.chdir('/home/rjackson/cmdv-rrm-anl/code/')
    import time_procedures
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import os
    from datetime import timedelta
    from copy import deepcopy
    import multidop

    # Get a Radar object given a time period in the CPOL dataset
    data_path_cpol = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_rapic/'
    out_file_path = '/home/rjackson/quicklook_plots/'
    out_data_path = '/home/rjackson/data/radar/grids/'

    # CPOL in lassen or rapic?
    cpol_format = 1    # 0 = lassen, 1 = rapic

    year_str = "%04d" % rad_time.year
    month_str = "%02d" % rad_time.month
    day_str = "%02d" % rad_time.day
    hour_str = "%02d" % rad_time.hour
    minute_str = "%02d" % rad_time.minute
    second_str = "%02d" % rad_time.second

    # Check to see if Cf/Radial file already exists...
    Grid = time_procedures.get_grid_from_cpol(rad_time)
    # Get reflectivity field
    reflectivity = Grid.fields['DT']['data']
    grid_shape = reflectivity.shape
 
    # Prepare shapes for statistical coverage product
    SCP0 = np.ma.zeros(grid_shape[0],)
    SCP10 = np.ma.zeros(grid_shape[0],)
    SCP20 = np.ma.zeros(grid_shape[0],) 
    SCP30 = np.ma.zeros(grid_shape[0],)
    SCP40 = np.ma.zeros(grid_shape[0],)
    total_points = grid_shape[1]*grid_shape[2]

    # Derive statistical coverage product
    for levels in range(0,grid_shape[0]):
        numgt0 = len(np.where(reflectivity[levels] > 0))
        SCP0[levels] = numgt0/total_points*100
        numgt10 = len(np.where(reflectivity[levels] > 0))
        SCP10[levels] = numgt10/total_points*100
        numgt20 = len(np.where(reflectivity[levels] > 0))
        SCP20[levels] = numgt20/total_points*100     
        numgt30 = len(np.where(reflectivity[levels] > 0))
        SCP30[levels] = numgt30/total_points*100    
        numgt40 = len(np.where(reflectivity[levels] > 0))
        SCP40[levels] = numgt40/total_points*100

    
    return SCP0, SCP10, SCP20, SCP30, SCP40
  
times = time_procedures.get_grid_times_cpol(start_year, 
                                            start_month,
                                            start_day,
                                            start_hour, 
                                            start_minute,
                                            end_year, 
                                            end_month,
                                            end_day,
                                            end_hour, 
                                            end_minute,
                                           )


# Go through all of the scans
#for rad_time in times:
#    display_time(rad_time)

# Get iPython cluster
#state = 0
#while state == 0:
#    try:
#        My_Cluster = Client()
#        My_View = My_Cluster[:]
#        state = 1
#    except:
#        state = 0
#        print('Cluster not ready for me')
#        sleep(10)

#Turn off blocking so all engines can work async
#My_View.block = False

#on all engines do an import of Py-ART
#My_View.execute('import matplotlib')
#My_View.execute('matplotlib.use("agg")')
#My_View.execute('import pyart')
#My_View.execute('import numpy as np')
#My_View.execute('from time_procedures import get_radar_from_cpol_rapic')
#t1 = time.time()
num_levels = 40
vert_levels = np.arange(0.5, 20.5, 0.5)
SCP0 = np.zeros((len(times), num_levels))
SCP10 = np.zeros((len(times), num_levels))
SCP20 = np.zeros((len(times), num_levels))
SCP30 = np.zeros((len(times), num_levels))
SCP40 = np.zeros((len(times), num_levels))
i = 0
for timer in times:
    SCP0[i,:], SCP10[i,:], SCP20[i,:], SCP30[i,:], SCP40[i,:] = SCP(timer)
    i = i + 1

# Output results into netCDF file
out_netcdf = Dataset('SCP.cdf', mode='w')
out_netcdf.createDimension('time', len(times))
out_netcdf.createDimension('levels', num_levels)

levels_file = out_netcdf.createVariable('levels', 'f8', ('levels',))
levels_file.long_name = 'Height in km'
levels_file.units = 'km'
levels_file[:] = vert_levels

SCP0_file = out_netcdf.createVariable('SCP0', 'f8', ('time','levels'))
SCP0_file.long_name = '% of level > 0 dBZ'
SCP0_file.units = '%'
SCP0_file = SCP0

SCP10_file = out_netcdf.createVariable('SCP10', 'f8', ('time','levels'))
SCP10_file.long_name = '% of level > 10 dBZ'
SCP10_file.units = '%'
SCP10_file = SCP10

SCP20_file = out_netcdf.createVariable('SCP20', 'f8', ('time','levels'))
SCP20_file.long_name = '% of level > 20 dBZ'
SCP20_file.units = '%'
SCP20_file = SCP20

SCP30_file = out_netcdf.createVariable('SCP30', 'f8', ('time','levels'))
SCP30_file.long_name = '% of level > 30 dBZ'
SCP30_file.units = '%'
SCP30_file = SCP30

SCP40_file = out_netcdf.createVariable('SCP40', 'f8', ('time','levels'))
SCP40_file.long_name = '% of level > 40 dBZ'
SCP40_file.units = '%'
SCP40_file = SCP40

years = out_netcdf.createVariable('years', 'i4', ('time',))
years.long_name = 'Year'
years.units = 'YYYY'
years = timer.year

days = out_netcdf.createVariable('days', 'i4', ('time',))
days.long_name = 'Day'
days.units = 'DD'
days = timer.day

months = out_netcdf.createVariable('months', 'i4', ('time',))
months.long_name = 'Month'
months.units = 'MM'
months = timer.month

hours = out_netcdf.createVariable('hours', 'i4', ('time',))
hours.long_name = 'Hour'
hours.units = 'HR'
hours = timer.hour

minutes = out_netcdf.createVariable('minutes', 'i4', ('time',))
minutes.long_name = 'Minute'
minutes.units = 'MM'
minutes = timer.minute

out_netcdf.close()

#Map the code and input to all workers
#result = My_View.map_async(SCP, times)

#Reduce the result to get a list of output
#qvps = result.get()
#tt = time.time() - t1
#print(tt)
#print(tt/len(times))


    
