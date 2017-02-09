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

# File paths
berr_data_file_path = '/lcrc/group/earthscience/radar/stage/radar_disk_two/berr_rapic/'
data_path_cpol = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_rapic/'
out_file_path = '/lcrc/group/earthscience/rjackson/quicklook_plots_no_dealiasing/'

## Berrima - 2009-2011 (new format), 2005-2005 (old format)
start_year = 2005
start_month = 1
start_day = 23
start_hour = 1
start_minute = 1

end_year = 2007
end_month = 5
end_day = 24
end_hour = 1
end_minute = 2

def display_time(rad_date):
    import pyart
    import matplotlib
    import sys
    sys.path.append('/home/rjackson/cmdv-rrm-anl/code/')
    import time_procedures
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import os
    from datetime import timedelta

    # Get a Radar object given a time period in the CPOL dataset
    data_path_cpol = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_rapic/'
    out_file_path = '/lcrc/group/earthscience/rjackson/quicklook_plots_no_dealias/'
    out_data_path = '/lcrc/group/earthscience/rjackson/cpol/'

    # CPOL in lassen or rapic?
    cpol_format = 1    # 0 = lassen, 1 = rapic
    
    def get_radar_from_cpol_rapic(time):
        from datetime import timedelta, datetime
        year_str = "%04d" % time.year
        month_str = "%02d" % time.month
        day_str = "%02d" % time.day
        hour_str = "%02d" % time.hour
        minute_str = "%02d" % time.minute
        second_str = "%02d" % time.second
        if(time.year == 2009 or (time.year == 2010 and
                                 time.month < 6)):
            dir_str = 'cpol_0910/rapic/'
        else:
            dir_str = 'cpol_1011/rapic/'

        file_name_str = (data_path_cpol +
                         dir_str + 
                         year_str +
                         month_str +
                         day_str +
                         hour_str +
                         minute_str +
                        'Gunn_Pt' +
                        '.rapic')
        radar = pyart.aux_io.read_radx(file_name_str)
        return radar 
    
    one_day_later = rad_date+timedelta(days=1)
    times, dates = time_procedures.get_radar_times_cpol_cfradial(rad_date.year, 
                                                                 rad_date.month,
                                                                 rad_date.day,
                                                                 1, 
                                                                 1,
                                                                 one_day_later.year, 
                                                                 one_day_later.month,
                                                                 one_day_later.day,
                                                                 0, 
                                                                 1,
                                                                 )
    print(times)
    for rad_time in times:
        year_str = "%04d" % rad_time.year
        month_str = "%02d" % rad_time.month
        day_str = "%02d" % rad_time.day
        hour_str = "%02d" % rad_time.hour
        minute_str = "%02d" % rad_time.minute
        second_str = "%02d" % rad_time.second

        # Check to see if Cf/Radial file already exists...
        out_path = (out_data_path +
	           '/' +
                    year_str +
                    '/' +
	            month_str +
                    '/' +
                    day_str +
	            '/')
        if not os.path.exists(out_path):
            os.makedirs(out_path)

        try:  
            radar = time_procedures.get_radar_from_cpol_cfradial(rad_time)
            
            if(rad_time.year > 2008):
                ref_field = 'Refl'
                vel_field = 'Vel'
            else:
                ref_field = 'reflectivity'
                vel_field = 'velocity'
   
            
            # Dealias velocities
            #gatefilter = pyart.correct.despeckle.despeckle_field(radar,
            #                                                     vel_field)
            #gatefilter.exclude_below(ref_field, 0)
                 
                
    
            
            # Save to Cf/Radial file
            
            out_path = (out_file_path +
                        '/' +
                        year_str +
                        '/' +
                        month_str +
                        '/' +
                        day_str +
                        '/')
	    
            if not os.path.exists(out_path):
                try:
                    os.makedirs(out_path)
                except: 
                    print('Not making directories')

            out_file = hour_str + minute_str + '.png'
            plt.figure(figsize=(7,14))
            plt.subplot(211)
            display = pyart.graph.RadarMapDisplay(radar)
            display.plot_ppi(ref_field, 
                             sweep=0, 
                             cmap=pyart.graph.cm.NWSRef,
                             vmin=0, 
                             vmax=70)
            plt.subplot(212)
            display = pyart.graph.RadarMapDisplay(radar)
            display.plot_ppi(vel_field, 
                             sweep=0,
                             cmap=pyart.graph.cm.NWSVel,
                             vmin=-30,
                             vmax=30)
            plt.savefig(out_path + out_file)
            plt.close() 
        except:
            print('Skipping corrupt time' +
                  year_str + 
                  '-' +
                  month_str + 
                  ' ' + 
                  hour_str + 
                  ':' +
                  minute_str)
            import sys
            print "'Error: '",  sys.exc_info()[0], sys.exc_info()[1]
times,dates = time_procedures.get_radar_times_cpol_cfradial(start_year, 
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
print(dates)

serial = 0
if(serial == 0):
    # Go through all of the scans
    #for rad_time in times:
    #    display_time(rad_time)

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

    #Map the code and input to all workers
    result = My_View.map_async(display_time, dates)

    #Reduce the result to get a list of output
    qvps = result.get()
    tt = time.time() - t1
    print(tt)
    print(tt/len(times))
else:
    for date in dates:
        display_time(date)


    
