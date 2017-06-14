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
import sys

# File paths
berr_data_file_path = '/lcrc/group/earthscience/radar/stage/radar_disk_two/berr_rapic/'
data_path_cpol = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_rapic/'
out_file_path = '/lcrc/group/earthscience/rjackson/quicklook_plots/'

## Berrima - 2009-2011 (new format), 2005-2005 (old format)
start_year = 2006
start_month = 1
start_day = 1
start_hour = 0
start_minute = 1

end_year = 2006
end_month = 1
end_day = 1
end_hour = 4
end_minute = 2

def grid_time(rad_time):
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
    import time

    # Get a Radar object given a time period in the CPOL dataset
    data_path_cpol = '/lcrc/group/earthscience/radar/stage/radar_disk_two/cpol_rapic/'
    out_file_path = '/lcrc/group/earthscience/rjackson/quicklook_plots/'
    out_data_path = '/lcrc/group/earthscience/rjackson/data/radar/grids/'

    # CPOL in lassen or rapic?
    cpol_format = 1    # 0 = lassen, 1 = rapic

    year_str = "%04d" % rad_time.year
    month_str = "%02d" % rad_time.month
    day_str = "%02d" % rad_time.day
    hour_str = "%02d" % rad_time.hour
    minute_str = "%02d" % rad_time.minute
    second_str = "%02d" % rad_time.second

    # Check to see if Cf/Radial file already exists...
    out_path = (out_data_path +
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
            print(out_path + ' directory already exists!')

    cpol_grid_name = (out_path + 'cpol_' 
                               + year_str 
                               + month_str
                               + day_str
                               + hour_str
                               + minute_str + '.nc')
    print(cpol_grid_name)
    try:
            Radar = time_procedures.get_radar_from_cpol_cfradial(rad_time)
            # Skip single sweep scans
            if(Radar.nsweeps == 1):
                return
            if(rad_time.year > 2007):
                ref_field = 'Refl'
                vel_field = 'Vel'
            else:
                ref_field = 'reflectivity'
                vel_field = 'velocity'

            print('Calculating texture....')
            
        	    
            # Apply gatefilter based on velocity and despeckling
            bt = time.time()
            gatefilter_Gunn = pyart.correct.despeckle_field(Radar, 
                                                            ref_field, 
                                                            size=6)
            nyq_Gunn = Radar.instrument_parameters['nyquist_velocity']['data'][0]
            try:
                data = ndimage.filters.generic_filter(Radar.fields['Vel']['data'],
                                                      pyart.util.interval_std, size = (3,3),
                                                      extra_arguments = (-nyq_Gunn, nyq_Gunn))
                filtered_data = ndimage.filters.median_filter(data, size = (3,3))
                texture_field = pyart.config.get_metadata('Vel')
                texture_field['data'] = filtered_data
                Radar.add_field('velocity_texture', texture_field, replace_existing = True)
            
                gatefilter = pyart.correct.GateFilter(Radar)
                #gatefilter.exclude_below('Refl', -5)
                #gatefilter.exclude_below('velocity_texture', 3)
                gatefilter = pyart.correct.despeckle.despeckle_field(Radar,
                                                                     ref_field,
                                                                     size=6,
                                                                     gatefilter=gatefilter)
            except:
                try:
                    data = ndimage.filters.generic_filter(Radar.fields['velocity']['data'],
                                                          pyart.util.interval_std, size = (3,3),
                                                          extra_arguments = (-nyq_Gunn, nyq_Gunn))
                    filtered_data = ndimage.filters.median_filter(data, size = (3,3))
                    texture_field = pyart.config.get_metadata('velocity')
                    texture_field['data'] = filtered_data
                    Radar.add_field('velocity_texture', texture_field, replace_existing = True)
            
                    gatefilter = pyart.correct.GateFilter(Radar)
                    #gatefilter.exclude_below('reflectivity', -5)
                    #gatefilter.exclude_below('velocity_texture', 3)
                    gatefilter = pyart.correct.despeckle.despeckle_field(Radar,
                                                                         ref_field,
                                                                         size=6,
                                                                         gatefilter=gatefilter)
                except:
                    print('Unrecognized reflectivity/velocity field names! Skipping...')
                    return
            print(str(time.time() - bt) + ' seconds to filter')
            # Change variable names to DT (reflectivity) and VT (velocity) expected by multidop
            # If you needed to dealias or perform other corrections,
            # this would be the time to start doing that.
            # Both datasets already have aliasing corrections
            cp = deepcopy(Radar.fields[ref_field]['data'])
            #texture = Radar.fields['velocity_texture']['data']
            Radar.add_field_like(ref_field, 'DT', cp, replace_existing=True)
            #cp = deepcopy(Radar.fields['corrected_velocity']['data'])
            #Radar.add_field_like('corrected_velocity', 'VT', cp, replace_existing=True)
              
            # The analysis engine currently expects the "missing_value" attribute
            Radar.fields['DT']['missing_value'] = 1.0 * Radar.fields['DT']['_FillValue']
            #Radar.fields['VT']['missing_value'] = 1.0 * Radar.fields['VT']['_FillValue']

            # Grid the data to a Cartesian grid. The Dual doppler domain does not extend ~60 km 
            # from both radars, so no need to store more data than that. 
            grid_cpol = time_procedures.grid_radar(Radar, 
                                                   origin=(Radar.latitude['data'][0], 
                                                           Radar.longitude['data'][0]),
                                                   xlim=(-60000, 60000), 
                                                   ylim=(-60000, 60000), 
                                                   fields=['DT', 'velocity_texture'], 
                                                   min_radius=750.0, 
                                                   bsp=1.0, nb=1.5,
                                                   h_factor=3.0, 
                                                   gatefilter=gatefilter_Gunn,
                                                   zlim=(500, 20000), 
                                                   grid_shape=(40, 81, 111))
             
            

            # Save the input grids for later.
            pyart.io.write_grid(cpol_grid_name, grid_cpol)              
            
    except:
            import sys
            exc_type, exc_obj, exc_tb = sys.exc_info()
            print('Skipping corrupt time' +
                  year_str + 
                  '-' +
                  month_str + 
                  ' ' + 
                  hour_str + 
                  ':' +
                  minute_str)
            print('Exception: ' + 
                  str(sys.exc_info()[0]) + 
                  str(sys.exc_info()[1]) +
                  str(exc_tb.tb_lineno))

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

print(times)
# Go through all of the scans
#for rad_time in times:
#    grid_time(rad_time)

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
My_View.execute('from time_procedures import get_radar_from_cpol_rapic')
t1 = time.time()

#for timer in times:
#    grid_time(timer)

#Map the code and input to all workers
result = My_View.map_async(grid_time, times)

#Reduce the result to get a list of output
qvps=result.get()
print(result.stdout)
tt = time.time() - t1
print(tt)
print(tt/len(times))


    
