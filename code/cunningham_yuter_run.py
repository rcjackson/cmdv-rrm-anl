from cunningham_yuter_core import cunningham_yuter_grids
from ipyparallel import Client
import time_procedures
from time import sleep
import time
import numpy as np
from netCDF4 import Dataset

# Get the radars from given time.
# Input the range of dates and time wanted for the collection of images
start_year = 2006
start_day = 15
start_month = 1
start_hour = 1
start_minute = 1
start_second = 0

end_year = 2006
end_month = 2
end_day = 28
end_hour = 1
end_minute = 0
end_second = 0


# Make the sounding file for the input
times = time_procedures.get_radar_times_cpol(start_year, start_month, 
                                             start_day, start_hour, 
                                             start_minute, end_year,
                                             end_month, end_day, end_hour, 
                                             end_minute, minute_interval=0)

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
#   cunningham_yuter_grids(timer)

#Do Steiner convective classification (map to workers)
result = My_View.map_async(cunningham_yuter_grids, times)


#Reduce the result to get a list of output
qvps = result.get()
tt = time.time() - t1
print(tt)
print(tt/len(times))

# Load all of the grids
i = 0

# Number of scans to include in classification
nscans = 18 

for t in times:
    year_str = "%04d" % t.year
    day_str = "%02d" % t.day
    month_str = "%02d" % t.month 
    hour_str = "%02d" % t.hour
    minute_str = "%02d" % t.minute  
    fname = ('/home/rjackson/conv_stratiform/cpol_conv_strat' + year_str 
                                                              + month_str
                                                              + day_str
                                                              + hour_str
                                                              + minute_str +'.nc')
    Convective_netcdf = Dataset(fname, mode='r')

    # Load convective classification
    conv_strat = Convective_netcdf.variables['strat_conv'][:]
    try:
        precip[i,:,:] = conv_strat
    except:
        grid_size = conv_strat.shape
        precip = np.zeros((len(times), grid_size[0], grid_size[1])) 
        precip[i,:,:] = conv_strat
    i = i + 1

# Do classification for each nscans period surrounding the scan
precip_array_shape = precip.shape
organization = np.ma.zeros(precip.shape[0],)
scan_years = np.ma.zeros(precip.shape[0],)
scan_months = np.ma.zeros(precip.shape[0],)
scan_days = np.ma.zeros(precip.shape[0],)
scan_hours = np.ma.zeros(precip.shape[0],)
scan_minutes = np.ma.zeros(precip.shape[0],)

for i in range(0, precip_array_shape[0]):
    if(i - nscans/2 < 0 or i + nscans/2 > precip_array_shape[0]):
        organization[i] = float('nan')
    else:
        total_frequency = np.zeros((precip_array_shape[1], 
                                   precip_array_shape[2]))
        convective_frequency = np.zeros((precip_array_shape[1], 
                                        precip_array_shape[2]))
        total_changes = np.zeros((precip_array_shape[1], 
                                 precip_array_shape[2]))
        intermittency = np.zeros((precip_array_shape[1], 
                                 precip_array_shape[2]))
        modes = np.zeros((precip_array_shape[1], 
                         precip_array_shape[2]))
        precip_types = precip[int(i-nscans/2):int(i+nscans/2),:,:] 
   
        opportunities = nscans
        # Do Cunningham and Yuter for scans
        for j in range(0, precip_array_shape[1]):
            for k in range(0, precip_array_shape[2]):
               ptypes = np.squeeze(precip_types[:,j,k])
               total_frequency[j,k] = sum(np.where(ptypes > 0, 1, 0))
               if(total_frequency[j,k] > 0):
                   convective_frequency[j,k] = float(sum(np.where(ptypes == 2, 1, 0))) / float(total_frequency[j,k])
               changes = np.zeros(ptypes.shape)
               changes[1:] = np.diff(ptypes)
               changes[ptypes == 0] = 0
               total_changes[j,k] = sum(changes[ptypes > 0])
               intermittency[j,k] = float(total_changes[j,k])/float(opportunities)
               if(convective_frequency[j,k] < 0.33):
                   modes[j,k] = 1
               elif(convective_frequency[j,k] > 0.33 and
                    convective_frequency[j,k] < 0.66 and
                   intermittency[j,j] > 0.33):
                   modes[j,k] = 2
               elif(convective_frequency[j,k] > 0.66):
                   modes[j,k] = 3
               else:
                   modes[j,k] = 4
        
               # Little precipitation periods
               if(float(total_frequency[j,k])/float(opportunities) < 0.3):
                   modes[j,k] = 0
        
        # Output classification for period into variable
        total_points = len(np.where(modes > 0)[1])
        total_stratiform = len(np.where(modes == 1)[1])
        total_embedded = len(np.where(modes == 2)[1])
        total_convective = len(np.where(modes == 3)[1])
        total_other = len(np.where(modes == 4)[1])

        if(float(total_stratiform)/float(total_points) > 0.25):
            organization[i] = 0
        elif(float(total_convective)/float(total_points) > 0.25):
            organization[i] = 1
        elif(float(total_embedded)/float(total_points) > 0.25):    
            organization[i] = 2
        else:
            organization[i] = 3
  
    scan_years[i] = times[i].year
    scan_months[i] = times[i].month
    scan_days[i] = times[i].day
    scan_hours[i] = times[i].hour
    scan_minutes[i] = times[i].minute

# Output into netCDF file
out_netcdf = Dataset('CYClassification.cdf', mode='w')
out_netcdf.createDimension('time', len(organization))

org_file = out_netcdf.createVariable('organization', 'f4', ('time',))
org_file.long_name = 'Convective organization following Cunningham and Yuter [2013]'
org_file.units = '0 = stratiform, 1 = convective, 2 = embedded convective, 3 = other'
org_file[:] = organization

year_file = out_netcdf.createVariable('year', 'f4', ('time',))
year_file.long_name = 'Scan year'
year_file.units = 'YYYY'
year_file[:] = scan_years

month_file = out_netcdf.createVariable('month', 'f4', ('time',))
month_file.long_name = 'Scan month'
month_file.units = 'MM'
month_file[:] = scan_months

day_file = out_netcdf.createVariable('day', 'f4', ('time',))
day_file.long_name = 'Scan day'
day_file.units = 'MM'
day_file[:] = scan_days

hour_file = out_netcdf.createVariable('hour', 'f4', ('time',))
hour_file.long_name = 'Scan hour'
hour_file.units = 'HH'
hour_file[:] = scan_hours

minute_file = out_netcdf.createVariable('minute', 'f4', ('time',))
minute_file.long_name = 'Scan hour'
minute_file.units = 'HH'
minute_file[:] = scan_minutes

    
