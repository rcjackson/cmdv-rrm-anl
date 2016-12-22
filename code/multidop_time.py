import matplotlib
# Needed for Blues
matplotlib.use('agg')
from matplotlib import rcParams
from matplotlib import pyplot as plt
# PyART
import pyart
import gzip
import sys
from scipy import ndimage
import shutil, os
from datetime import timedelta, datetime
import numpy as np
import tempfile
import re
import time
from copy import deepcopy
from IPython.display import Image, display
import multidop
import math
from ipyparallel import Client
from time import sleep
import time_procedures

def do_multidop_for_time(frame_time):
    import matplotlib
    # Needed for Blues
    matplotlib.use('agg')
    from matplotlib import rcParams
    
    # PyART
    import sys
    from scipy import ndimage
    import shutil, os
    from datetime import timedelta, datetime
    import tempfile
    import glob
    import re
    import time
    from copy import deepcopy
    from IPython.display import Image, display
    import multidop
    import math
    from netCDF4 import Dataset
    from datetime import timedelta, datetime
    import time_procedures
    import os
    import numpy as np

    year_str = "%04d" % frame_time.year
    day_str = "%02d" % frame_time.day
    month_str = "%02d" % frame_time.month 
    hour_str = "%02d" % frame_time.hour
    minute_str = "%02d" % frame_time.minute  
    fname = (time_procedures.out_data_path + 'ddop/cf_compliant_grid' 
                                           + year_str 
                                           + month_str
                                           + day_str
                                           + hour_str
                                           + minute_str +'.nc')
    print('Does file ' + fname + ' exist?')
    if(not os.path.isfile(fname)):
        one_day_ago = frame_time-timedelta(days=1, minutes=1)
        sounding_times = time_procedures.get_sounding_times(one_day_ago.year,
                                                            one_day_ago.month,
                                                            one_day_ago.day,
                                                            one_day_ago.hour,
                                                            one_day_ago.minute,
                                                            frame_time.year,
                                                            frame_time.month,
                                                            frame_time.day,
                                                            frame_time.hour,
                                                            frame_time.minute,
                                                            minute_interval=60)
        print(sounding_times)
        sounding_time = sounding_times[len(sounding_times)-1]
        Sounding_netcdf = time_procedures.get_sounding(sounding_time)

        # Convert timestamps to datetime format
        Time = Sounding_netcdf.variables['time_offset'][:]
        base_time = Sounding_netcdf.variables['base_time'][:]
        alt = Sounding_netcdf.variables['alt'][:]
        u = Sounding_netcdf.variables['u_wind'][:]
        v = Sounding_netcdf.variables['v_wind'][:]
        base_timestamp = timedelta(seconds=float(base_time)) + datetime(1970,1,1)
        Sounding_netcdf.close()
    	
        time_delta = datetime(frame_time.year,
                              frame_time.month,
                              frame_time.day,
                              frame_time.hour,
                              frame_time.minute,
                              frame_time.second) - base_timestamp
        seconds_in_file = time_delta.days*(24*60*60) + time_delta.seconds
        
        one_minute_later = frame_time+timedelta(minutes=1)
        ten_minutes_ago = frame_time-timedelta(minutes=10)
        one_minute_earlier = frame_time-timedelta(minutes=1)

        times_berr = time_procedures.get_radar_times_berr(one_minute_earlier.year,
                                                          one_minute_earlier.month,
                                                          one_minute_earlier.day,
                                                          one_minute_earlier.hour,
                                                          one_minute_earlier.minute,                                       
                                                          one_minute_later.year,
                                                          one_minute_later.month,
                                                          one_minute_later.day,
                                                          one_minute_later.hour,
                                                          one_minute_later.minute,
                                                          minute_interval=0)
        print(times_berr)
            
        sounding_file_name = (time_procedures.out_data_path + 'soundings/'
                                                            + year_str 
                                                            + month_str
                                                            + day_str
                                                            + hour_str
                                                            + minute_str)
                                              

        file = open(sounding_file_name, 'w')
 
        # Take 1000 evenly spaced levels from the sounding and place them into the file
        us = u[u > -75]
        vs = v[u > -75]
        alts = alt[u > -75]
        step = int(math.floor(len(us)/500))

        for i in range(0, len(us), step):
            input_string = (str(alts[i]) + ' ' + str(us[i]) + ' ' + str(vs[i]) + '\n')
            file.write(input_string)

        file.close()
        # Calculate texture of velocity field for Berrima and CPOL
        # if previous frame is not available, just use (u,v) = 0
    

        Radar = time_procedures.get_radar_from_cpol(frame_time)
        try:
            Radar_berr = time_procedures.get_radar_from_berr(times_berr[0])
        except:
            print('Cannot find matching time from Berrima radar, skipping')
            return
	
	bt = time.time()
        print('Calculating texture....')
        nyq_Gunn = Radar.instrument_parameters['nyquist_velocity']['data'][0]
        nyq_Berr = Radar_berr.instrument_parameters['nyquist_velocity']['data'][0]
        data = ndimage.filters.generic_filter(Radar.fields['velocity']['data'],
                                             pyart.util.interval_std, size = (4,4),
                                             extra_arguments = (-nyq_Gunn, nyq_Gunn))
        filtered_data = ndimage.filters.median_filter(data, size = (4,4))
        texture_field = pyart.config.get_metadata('velocity')
        texture_field['data'] = filtered_data
        Radar.add_field('velocity_texture', texture_field, replace_existing = True)
        data = ndimage.filters.generic_filter(Radar_berr.fields['velocity']['data'],
                                              pyart.util.interval_std, size = (4,4),
                                              extra_arguments = (-nyq_Gunn, nyq_Gunn))
        filtered_data = ndimage.filters.median_filter(data, size = (4,4))
        texture_field = pyart.config.get_metadata('velocity')
        texture_field['data'] = filtered_data
        Radar_berr.add_field('velocity_texture', texture_field, replace_existing = True)
        print('Done!')
	print((time.time()-bt)/60.0, 'minutes to process')
	    
        # Apply gatefilter based on velocity and despeckling
        gatefilter_Gunn = pyart.correct.despeckle_field(Radar, 
                                                        'corrected_reflectivity', 
                                                        size=6)
        gatefilter_Gunn.exclude_above('velocity_texture', 4)
	gatefilter_Gunn.exclude_below('corrected_reflectivity', 1)

        gatefilter_Berr = pyart.correct.despeckle_field(Radar_berr, 
                                                        'corrected_reflectivity', 
                                                        size=6)
        gatefilter_Berr.exclude_above('velocity_texture', 4)
	gatefilter_Gunn.exclude_below('corrected_reflectivity', 1)

        # Change variable names to DT (reflectivity) and VT (velocity) expected by multidop
        # If you needed to dealias or perform other corrections,
        # this would be the time to start doing that.
        # Both datasets already have aliasing corrections
        cp = deepcopy(Radar.fields['corrected_reflectivity']['data'])
        texture = Radar.fields['velocity_texture']['data']
        Radar.add_field_like('corrected_reflectivity', 'DT', cp, replace_existing=True)
        cp = deepcopy(Radar.fields['corrected_velocity']['data'])
        Radar.add_field_like('corrected_velocity', 'VT', cp, replace_existing=True)
       
        cp = deepcopy(Radar_berr.fields['corrected_reflectivity']['data'])
        Radar_berr.add_field_like('corrected_reflectivity', 'DT', 
                                  cp, replace_existing=True)
        cp = deepcopy(Radar_berr.fields['corrected_velocity']['data'])
        Radar_berr.add_field_like('corrected_velocity', 'VT', 
                                  cp, replace_existing=True)
        
        # The analysis engine currently expects the "missing_value" attribute
        Radar.fields['DT']['missing_value'] = 1.0 * Radar.fields['DT']['_FillValue']
        Radar_berr.fields['DT']['missing_value'] = 1.0 * Radar_berr.fields['DT']['_FillValue']
        Radar.fields['VT']['missing_value'] = 1.0 * Radar.fields['VT']['_FillValue']
        Radar_berr.fields['VT']['missing_value'] = 1.0 * Radar_berr.fields['VT']['_FillValue']

        # Grid the data to a Cartesian grid. The Dual doppler domain does not extend ~60 km 
        # from both radars, so no need to store more data than that. 
        grid_cpol = time_procedures.grid_radar(Radar, 
                                               origin=(Radar.latitude['data'][0], 
                                               Radar.longitude['data'][0]),
                                               xlim=(-60000, 50000), 
                                               ylim=(-50000, 30000), 
                                               fields=['DT', 'VT'], 
                                               min_radius=750.0, 
                                               bsp=1.0, nb=1.5,
                                               h_factor=3.0, 
                                               gatefilter=gatefilter_Gunn,
                                               zlim=(500, 16000), 
                                               grid_shape=(32, 81, 111))
        grid_Berr = time_procedures.grid_radar(Radar_berr, 
                                               origin=(Radar.latitude['data'][0], 
                                               Radar.longitude['data'][0]),
                                               fields=['DT', 'VT'],
                                               xlim=(-60000, 50000), 
                                               ylim=(-50000, 30000), 
                                               zlim=(500, 16000), 
                                               min_radius=750.0,  
                                               grid_shape=(32, 81, 111), 
                                               gatefilter=gatefilter_Berr,
                                               bsp=1.0, nb=1.5, 
                                               h_factor=4.0)

	grid_Berr.fields['DT']['data'] = grid_cpol.fields['DT']['data']
        # The analysis engine requires azimuth and elevation to be part of the grid.
        # This information is computed from the grid geometry.
        grid_cpol = multidop.angles.add_azimuth_as_field(grid_cpol)
        grid_Berr = multidop.angles.add_azimuth_as_field(grid_Berr)
        grid_cpol = multidop.angles.add_elevation_as_field(grid_cpol)
        grid_Berr = multidop.angles.add_elevation_as_field(grid_Berr)

        cpol_grid_name = (time_procedures.out_data_path + 'cpol/cpol_' 
                                                        + year_str 
                                                        + month_str
                                                        + day_str
                                                        + hour_str
                                                        + minute_str +'.nc')

        berr_grid_name = (time_procedures.out_data_path + 'berr/berr_' 
                                                        + year_str 
                                                        + month_str
                                                        + day_str
                                                        + hour_str
                                                        + minute_str +'.nc')
        # Save the input grids for later.
        pyart.io.write_grid(cpol_grid_name, grid_cpol)
        pyart.io.write_grid(berr_grid_name, grid_Berr)
  
        # Load previous time period for storm motion (use 0 if no previous frame)
        try:
            Radar_prev = get_radar_from_cpol(ten_minutes_ago)

            print('Calculating storm motion....')
            nyq_Gunn = Radar_prev.instrument_parameters['nyquist_velocity']['data'][0]
            data = ndimage.filters.generic_filter(Radar_prev.fields['velocity']['data'],
                                                  pyart.util.interval_std, size = (4,4),
                                                  extra_arguments = (-nyq_Gunn, nyq_Gunn))
            filtered_data = ndimage.filters.median_filter(data, size = (4,4))
            texture_field = pyart.config.get_metadata('velocity')
            texture_field['data'] = filtered_data

            print('Gridding previous frame...')
            cp = deepcopy(Radar_prev.fields['corrected_reflectivity']['data'])
            cp = np.ma.masked_where(texture_field['data'] > 4, cp)
            Radar_prev.add_field_like('corrected_reflectivity', 'DT', cp, 
                                      replace_existing=True)
            grid_prev = time_procedures.grid_radar(Radar_prev, 
                                                   origin=(Radar_prev.latitude['data'][0], 
                                                   Radar_prev.longitude['data'][0]),
                                                   xlim=(-60000, 60000), 
                                                   ylim=(-50000, 30000), 
                                                   fields=['DT'],
                                                   zlim=(500, 16000), 
                                                   grid_shape=(32, 121, 81))
            (vt,ut) = pyart.retrieve.grid_displacement_pc(grid_prev, grid_cpol, 
                                                          'DT', 9, 
                                                          return_value='velocity')
        except:
            (vt,ut) = (0,0)

        # You don't have to define everything. 
        # Most of these keywords are default values.
        # If you don't define something the program will provide a default value.
        # Check parameters.py for what keyword default values are.
        calc_file_name = (time_procedures.out_data_path + '/dda_files/cpol_calc' 
                                                        + year_str 
                                                        + month_str
                                                        + day_str
                                                        + hour_str
                                                        + minute_str +'.dda')
    
        frprmn_out_name = (time_procedures.out_data_path + '/dda_files/frprmn_out' 
                                                        + year_str 
                                                        + month_str
                                                        + day_str
                                                        + hour_str
                                                        + minute_str +'.nc')
        localfile = tempfile.NamedTemporaryFile()
        pd = {'dir': './',
              'x': [-60000.0, 1000.0, 111],   # start, step, max = min + (steps-1)
              'y': [-50000.0, 1000.0, 81],
              'z': [500.0, 500.0,  32],
              'grid': [grid_cpol.origin_longitude['data'][0], 
                       grid_cpol.origin_latitude['data'][0], 
                       50.0],
              'files': [berr_grid_name,
                      cpol_grid_name],
              'radar_names': ['Berrima', 'CPOL'],
              'refl': 'DT',  # Name of reflectivity field. Must be common between radars.
              'vt': 'VT',  # Name of velocity field. Must be common between radars.
              'bgfile': sounding_file_name, # Name of sounding file
              'writeout': localfile.name, # Name of output grid file
              'frprmn_out': frprmn_out_name,
              'min_cba': 30.0,  # Minimum beam-crossing angle
              'calc_params': calc_file_name, # .dda file for parameters 
                                             # related to minimization routine
              'anel': 1, # 0 = Boussinesq approximation  1 = anelastic 
              'laplace': 0, # 0 = 1st order derivatives for smoothing, 1 = second
              'read_dataweights': 2, # 0 = calculate data constraint weights/output, 
                                     # 1 = read from file, 2 = weigh all equally
              'max_dist': 10.0, # How much distance analysis and observational 
                                # grid must match in m
              'cutoff': 0.0, # Deny observations below this level from analysis (m)
              'UT': ut, # U of prescribed storm motion vector
              'VT': vt, # V of prescribed storm motion vector
              'output_error': 0, # 1 = output verification stats after each iteration
              'weak_height': -1, # Sounding height constraint weakened in regions 
                                 # > 10 dBZ below this height (-1 = disabled)
              'upper_bc': 1, # 1 = w = 0 as upper boundary condition, -1 = ignore
              'itmax_frprmn': [200, 10], # max iterations in frprmn function
              'itmax_dbrent': 200, # max iterations in dbrent function
              'C1b': 0.1,  # Data weighting factor
              'C2b': 1500.0,  # Mass continuity weighting factor
              'C3b': 0.0,  # Vorticity weighting factor
              'C4b': 100.0,  # Horizontal smoothing factor
              'C5b': 2.0,  # Vertical smoothing factor
              'C8b': 0.01,  # Sounding factor
              'vary_weights': 0,
              # Define filter with ONE of the following forms.
              # filter: none
              # filter: filter_frequency Leise nstep
              # filter: filter_frequency low-pass alpha
              'filter': ['60', 'Leise', '2'],
              # Coverage values for various combinations of radars.
              # Each line should provide the type of coverage value, radar count,
              # radar names, and the value, in the following form:
              #
              #   cvg_(""|opt|sub)_(bg|fil): integer radar1 radar2 ... boolean
              #
              # Radars are identified by the OPAWS/OBAN file name with grid data for that
              # radar. This must be just the base name, not the full path.
              #
              # For example:
              #
              #   cvg_opt_bg: SR1 SR2 1
              #
              # says that if SR1 SR2
              # both have data within max_dist meters of the point under consideration,
              # and an optimal beam crossing angle, then the point will receive a coverage
              # value of 1, i.e. point has coverage.
              #
              # "opt" means optimal beam crossing angle.
              # "sub" means suboptimal beam crossing angle.
              # "bg" means background coverage.
              # "fil" means filter coverage.
              # cvg_bg, cvg_fil, and sseq_trip do not require a radar count. (Beam crossing
              # angle is meaningless with one radar, so there is no opt or sub)
              #
              # If this file is being used, coverage values must be provided for all
              # combinations of radars.
              'cvg_opt_bg': [1, 1, 0],
              'cvg_sub_bg': [1, 1, 0],
              'cvg_opt_fil': [0, 0, 0],
              'cvg_sub_fil': [1, 1, 0],
              'cvg_bg': [1, 1, 0],
              'cvg_fil': [0, 0, 0],
              'sseq_trip': [1.0, 1.0, 0.0]
            }
        dda_file_name = (time_procedures.out_data_path + '/dda_files/cpol_test' 
                                                       + year_str 
                                                       + month_str
                                                       + day_str
                                                       + hour_str
                                                       + minute_str +'.dda')
  
        pf = multidop.parameters.ParamFile(pd, dda_file_name)
        pf = multidop.parameters.CalcParamFile(pd, calc_file_name)

        # Unfortunately, text output from the analysis engine (DDA) will not display
        # until after the program completes. Expect this step to take several minutes.
        bt = time.time()
        #multidop.execute.do_analysis(dda_file_name, cmd_path='/home/rjackson/multidop/src/DDA')
        multidop.execute.run_command('./DDA ' + dda_file_name)
        print((time.time()-bt)/60.0, 'minutes to process')

        # Baseline output is not CF or Py-ART compliant. This function fixes that.
        # This is why we wrote the original output to a tempfile that can be safely removed.
        # The final grid will have all wind solutions outside the coverage region masked.
     
        final_grid = multidop.grid_io.make_new_grid([grid_cpol, grid_Berr], localfile.name)
        final_grid.write(fname)
        localfile.close()
    else:
        print('DDA grid already exists...skipping.')

