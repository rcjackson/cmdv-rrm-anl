from multidop_time import do_multidop_for_time
from ipyparallel import Client
import time_procedures
from time import sleep
import time
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import sys
import os
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

start_time = datetime(start_year, start_month, start_day,
                      start_hour, start_minute, start_second)
end_time = datetime(end_year, end_month, end_day,
                    end_hour, end_minute, end_second)

scp_netcdf_file_path = 'SCP.cdf'

# To save core hours, do not process clear air periods
scp_netcdf = Dataset(scp_netcdf_file_path, mode='r')

years = scp_netcdf.variables['years'][:]
months = scp_netcdf.variables['months'][:]
days = scp_netcdf.variables['days'][:]
hours = scp_netcdf.variables['hours'][:]
minutes = scp_netcdf.variables['minutes'][:]
SCP20 = scp_netcdf.variables['SCP20'][:,:]
levels = scp_netcdf.variables['levels'][:]
maxz = scp_netcdf.variables['maxz'][:]
scp_netcdf.close()

times = []
for i in range(0,len(years)):
    temp_date = datetime(year=years[i],
                         month=months[i],
                         day=days[i],
                         hour=hours[i],
                         minute=minutes[i]-1,
                         )
    print(minutes[i])
    if(SCP20[i,2] > 0.1 and temp_date >= start_time and temp_date <= end_time):
        # Check for file existence
        year_str = "%04d" % temp_date.year
        day_str = "%02d" % temp_date.day
        month_str = "%02d" % temp_date.month 
        hour_str = "%02d" % temp_date.hour
        minute_str = "%02d" % temp_date.minute  
        fname = (time_procedures.out_data_path + 'ddop/cf_compliant_grid' 
                                               + year_str 
                                               + month_str
                                               + day_str
                                               + hour_str
                                               + minute_str +'.nc')    
        if(not os.path.isfile(fname)):
            times.append(temp_date)   

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
    result = My_View.map_async(do_multidop_for_time, times)


    #Reduce the result to get a list of output
    qvps = result.get()
    tt = time.time() - t1
    print(tt)
    print(tt/len(times))
else:
    for timer in times:
        do_multidop_for_time(timer)


    
