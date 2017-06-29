import matplotlib
from matplotlib import rcParams
from matplotlib import pyplot as plt
matplotlib.use('agg')
import pyart
import gzip
import sys
from scipy import ndimage
import shutil
import os
from datetime import timedelta, datetime
import numpy as np
import tempfile
import re
import time
from copy import deepcopy
import multidop
import math
import time_procedures
import sys


def do_multidop_for_time(frame_time):
    year_str = "%04d" % frame_time.year
    day_str = "%02d" % frame_time.day
    month_str = "%02d" % frame_time.month
    hour_str = "%02d" % frame_time.hour
    minute_str = "%02d" % frame_time.minute
    fname = (time_procedures.out_data_path + 'ddop/cf_compliant_grid' +
             year_str + month_str + day_str + hour_str +
             minute_str + '.nc')
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
        sounding_time = sounding_times[len(sounding_times)-1]
        Sounding_netcdf = time_procedures.get_sounding(sounding_time)

        # Convert timestamps to datetime format
        Time = Sounding_netcdf.variables['time_offset'][:]
        base_time = Sounding_netcdf.variables['base_time'][:]
        alt = Sounding_netcdf.variables['alt'][:]
        u = Sounding_netcdf.variables['u_wind'][:]
        v = Sounding_netcdf.variables['v_wind'][:]
        base_timestamp = timedelta(seconds=float(base_time)) + datetime(
                                                               1970, 1, 1)
        Sounding_netcdf.close()

        time_delta = datetime(frame_time.year,
                              frame_time.month,
                              frame_time.day,
                              frame_time.hour,
                              frame_time.minute,
                              frame_time.second) - base_timestamp
        seconds_in_file = time_delta.days*(24*60*60) + time_delta.seconds

        one_minute_later = frame_time+timedelta(minutes=5)
        ten_minutes_ago = frame_time-timedelta(minutes=10)
        one_minute_earlier = frame_time-timedelta(minutes=5)

        times_berr, dates = time_procedures.get_radar_times_berr(
                            one_minute_earlier.year, one_minute_earlier.month,
                            one_minute_earlier.day, one_minute_earlier.hour,
                            one_minute_earlier.minute, one_minute_later.year,
                            one_minute_later.month, one_minute_later.day,
                            one_minute_later.hour, one_minute_later.minute,
                            minute_interval=0)

        # Load sounding data (ARM format assumed)
        sounding_file_name = (time_procedures.out_data_path + 'soundings/' +
                              year_str + month_str + day_str + hour_str +
                              minute_str)
        file = open(sounding_file_name, 'w')

        """ Take 1000 evenly spaced levels from the sounding
             and place them into the file. Multidop needs a file
             that is space separated with each row being:
             altitude [m], u [m/s], w [m/s] """

        # Do not include invalid/missing entries
        us = u[u > -75]
        vs = v[u > -75]
        alts = alt[u > -75]
        step = int(math.floor(len(us)/500))
        if(step > 0):
            for i in range(0, len(us), step):
                input_string = (str(alts[i]) + ' ' + str(us[i]) +
                                ' ' + str(vs[i]) + '\n')
                file.write(input_string)

        # If baloon popped below 15 km (approximate tropopause),
        # then don't use sounding
        if(alts[-1] < 15000 or step == 0):
            use_sounding = 0
        else:
            use_sounding = 1

        file.close()
        # Calculate texture of velocity field for Berrima and CPOL
        # if previous frame is not available, just use (u,v) = 0
        try:
            Radar = time_procedures.get_radar_from_cpol(frame_time)
            if(Radar.nsweeps == 1):
                print('CPOL radar only has one sweep!')
                return
        except:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            print(('Could not load CPOL radar data file! ' +
                   str(sys.exc_info()[0]) + str(sys.exc_info()[1]) +
                   str(exc_tb.tb_lineno)))
            return

        try:
            Radar_berr = time_procedures.get_radar_from_berr(times_berr[0])
        except:
            print('Cannot find matching time from Berrima radar, skipping')
            return

        if(frame_time.year <= 2007):
            cpol_ref_field = 'reflectivity'
            cpol_vel_field = 'velocity'
        else:
            cpol_ref_field = 'Refl'
            cpol_vel_field = 'Vel'

        bt = time.time()
        print('Calculating texture....')
        try:
            texture_field = pyart.filters.calculate_velocity_texture(
                            Radar, vel_field=cpol_vel_field)
            Radar.add_field(
                'velocity_texture', texture_field, replace_existing=True)
            texture_field = pyart.filters.calculate_velocity_texture(
                            Radar_berr, vel_field='Vel')
            Radar_berr.add_field(
                'velocity_texture', texture_field, replace_existing=True)
            print('Done!')
            print((time.time()-bt)/60.0, 'minutes to process')
        except:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            print(('Could not calculate texture! ' +
                   str(sys.exc_info()[0]) + str(sys.exc_info()[1]) +
                   str(exc_tb.tb_lineno)))
            return

        # Apply gatefilter based on velocity and despeckling
        gatefilter_Gunn = pyart.correct.despeckle_field(Radar,
                                                        cpol_ref_field,
                                                        size=6)
        gatefilter_Gunn.exclude_above('velocity_texture', 3)
        gatefilter_Gunn.exclude_below(cpol_ref_field, 1)

        gatefilter_Berr = pyart.correct.despeckle_field(Radar_berr,
                                                        'Refl',
                                                        size=6)
        gatefilter_Berr.exclude_above('velocity_texture', 4)
        gatefilter_Berr.exclude_below('Refl', 1)

        # Change variable names to DT (reflectivity) and
        # VT (velocity) expected by multidop
        # This assumes datasets already have aliasing corrections
        cp = deepcopy(Radar.fields[cpol_ref_field]['data'])
        texture = Radar.fields['velocity_texture']['data']
        Radar.add_field_like(cpol_ref_field, 'DT', cp, replace_existing=True)
        try:
            cp = deepcopy(Radar.fields['corrected_velocity']['data'])
            Radar.add_field_like('corrected_velocity', 'VT',
                                 cp, replace_existing=True)
        except:
            print('No dealiased velocities from CPOL...skipping!')
            return

        cp = deepcopy(Radar_berr.fields['Refl']['data'])
        Radar_berr.add_field_like('Refl', 'DT',
                                  cp, replace_existing=True)
        try:
            cp = deepcopy(Radar_berr.fields['corrected_velocity']['data'])
            Radar_berr.add_field_like('corrected_velocity', 'VT',
                                      cp, replace_existing=True)
        except:
            print('No dealiased velocities from CPOL...skipping!')
            return
        # The analysis engine currently expects the "missing_value" attribute
        DT_cpol_FV = 1.0 * Radar.fields['DT']['_FillValue']
        DT_berr_FV = 1.0 * Radar_berr.fields['DT']['_FillValue']
        VT_cpol_FV = 1.0 * Radar.fields['VT']['_FillValue']
        VT_berr_FV = 1.0 * Radar_berr.fields['VT']['_FillValue']
        Radar.fields['DT']['missing_value'] = DT_cpol_FV
        Radar_berr.fields['DT']['missing_value'] = DT_berr_FV
        Radar.fields['VT']['missing_value'] = VT_cpol_FV
        Radar_berr.fields['VT']['missing_value'] = VT_berr_FV

        # Grid the data to a Cartesian grid. The Dual doppler domain
        # does not extend ~60 km from both radars, so no need to store
        # more data than that.
        origin = (Radar.latitude['data'][0], Radar.longitude['data'][0])
        grid_cpol = time_procedures.grid_radar(Radar,
                                               origin=origin,
                                               xlim=(-60000, 50000),
                                               ylim=(-50000, 30000),
                                               fields=['DT', 'VT'],
                                               min_radius=750.0,
                                               bsp=1.0, nb=1.5,
                                               h_factor=3.0,
                                               gatefilter=gatefilter_Gunn,
                                               zlim=(500, 20000),
                                               grid_shape=(40, 81, 111))
        grid_Berr = time_procedures.grid_radar(Radar_berr,
                                               origin=origin,
                                               fields=['DT', 'VT'],
                                               xlim=(-60000, 50000),
                                               ylim=(-50000, 30000),
                                               zlim=(500, 20000),
                                               min_radius=750.0,
                                               grid_shape=(40, 81, 111),
                                               gatefilter=gatefilter_Berr,
                                               bsp=1.0, nb=1.5,
                                               h_factor=3.0)

        # Berrima reflectivities are corrupt for many scans
        #  -- prefer CPOL reflectivities (remove for 2 good Z datasets)
        grid_Berr.fields['DT']['data'] = grid_cpol.fields['DT']['data']

        # The analysis engine requires azimuth and elevation
        # to be part of the grid. This information is computed from
        # the grid geometry.
        grid_cpol = multidop.angles.add_azimuth_as_field(grid_cpol)
        grid_Berr = multidop.angles.add_azimuth_as_field(grid_Berr)
        grid_cpol = multidop.angles.add_elevation_as_field(grid_cpol)
        grid_Berr = multidop.angles.add_elevation_as_field(grid_Berr)

        cpol_grid_name = (time_procedures.out_data_path + 'cpol/cpol_' +
                          year_str + month_str + day_str + hour_str +
                          minute_str + '.nc')

        berr_grid_name = (time_procedures.out_data_path + 'berr/berr_' +
                          year_str + month_str + day_str + hour_str +
                          minute_str + '.nc')

        # Save the input grids for later.
        pyart.io.write_grid(cpol_grid_name, grid_cpol)
        pyart.io.write_grid(berr_grid_name, grid_Berr)

        # Load previous time period for storm motion
        # (use 0 if no previous frame)
        try:
            Radar_prev = get_radar_from_cpol(ten_minutes_ago)

            print('Calculating storm motion....')
            texture_field = pyart.filters.calculate_velocity_texture(
                            Radar_prev, vel_field=cpol_vel_field)

            print('Gridding previous frame...')
            cp = deepcopy(Radar_prev.fields['corrected_reflectivity']['data'])
            cp = np.ma.masked_where(texture_field['data'] > 4, cp)
            Radar_prev.add_field_like('corrected_reflectivity', 'DT', cp,
                                      replace_existing=True)
            grid_prev = time_procedures.grid_radar(Radar_prev,
                                                   origin=origin,
                                                   xlim=(-60000, 60000),
                                                   ylim=(-50000, 30000),
                                                   fields=['DT'],
                                                   zlim=(500, 20000),
                                                   grid_shape=(40, 121, 81))
            (vt, ut) = pyart.retrieve.grid_displacement_pc(
                       grid_prev, grid_cpol, 'DT', 9,
                       return_value='corrected_velocity')
        except:
            (vt, ut) = (0, 0)

        # You don't have to define everything.
        # Most of these keywords are default values.
        # If you don't define something the program will
        # provide a default value.
        # Check parameters.py for what keyword default values are.
        calc_file_name = (time_procedures.out_data_path +
                          '/dda_files/cpol_calc' + year_str + month_str +
                          day_str + hour_str + minute_str + '.dda')

        frprmn_out_name = (time_procedures.out_data_path +
                           '/fprmn/frprmn_out' + year_str + month_str +
                           day_str + hour_str + minute_str + '.nc')
        localfile = tempfile.NamedTemporaryFile()

        # If sounding is available, favor constraint based on sounding
        # vs. data, otherwise favor data more
        if(use_sounding == 0):
            C8b = 0.0
            C1b = 1.0
            sounding_file_name = None
        else:
            C1b = 0.1
            C8b = 0.01

        pd = {'dir': './',
              'x': [-60000.0, 1000.0, 111],  # start, step, max=min+(steps-1)
              'y': [-50000.0, 1000.0, 81],
              'z': [500.0, 500.0,  40],
              'grid': [grid_cpol.origin_longitude['data'][0],
                       grid_cpol.origin_latitude['data'][0],
                       50.0],
              'files': [berr_grid_name, cpol_grid_name],
              'radar_names': ['Berrima', 'CPOL'],
              'refl': 'DT',  # Name of reflectivity field.
              'vt': 'VT',  # Name of velocity field.
              'bgfile': sounding_file_name,  # Name of sounding file
              'writeout': localfile.name,  # Name of output grid file
              'frprmn_out': frprmn_out_name,
              'min_cba': 30.0,  # Minimum beam-crossing angle
              'calc_params': calc_file_name,  # .dda file for parameters
                                              # related to minimization
              'anel': 1,  # 0 = Boussinesq approximation  1 = anelastic
              'laplace': 0,  # 0 = 1st order derivatives in smoothing, 1 = 2nd
              'read_dataweights': 2,  # 0 = calculate data constraint weights,
                                      # 1 = read from file, 2 = weigh equally
              'max_dist': 10.0,  # How much distance analysis and observational
                                 # grid must match in m
              'cutoff': 0.0,  # Deny observations below this level (m)
              'UT': ut,  # U of prescribed storm motion vector
              'VT': vt,  # V of prescribed storm motion vector
              'output_error': 0,  # 1 = output verification stats
              'weak_height': -1,  # Sounding height constraint weakened
                                  # > 10 dBZ below this height (-1 = off)
              'upper_bc': 1,  # 1 = w = 0 as upper boundary cond., -1 = ignore
              'itmax_frprmn': [200, 10],  # max iterations in frprmn function
              'itmax_dbrent': 200,  # max iterations in dbrent function
              'C1b': C1b,  # Data weighting factor
              'C2b': 1500.0,  # Mass continuity weighting factor
              'C3b': 0.0,  # Vorticity weighting factor
              'C4b': 75.0,  # Horizontal smoothing factor
              'C5b': 2.0,  # Vertical smoothing factor
              'C8b': C8b,  # Sounding factor
              'vary_weights': 0,
              # Define filter with ONE of the following forms.
              # filter: none
              # filter: filter_frequency Leise nstep
              # filter: filter_frequency low-pass alpha
              'filter': ['60', 'Leise', '2'],
              """ Coverage values for various combinations of radars.
              Each line should provide the type of coverage value, radar count,
              radar names, and the value, in the following form:

                  cvg_(""|opt|sub)_(bg|fil): integer radar1 radar2 ... boolean

              Radars are identified by the OPAWS/OBAN file name with grid
              data for that radar. This must be just the base name,
              not the full path.

              For example:

                 cvg_opt_bg: SR1 SR2 1

              says that if SR1 SR2 both have data within max_dist meters
              of the point under consideration, and an optimal beam crossing
              angle, then the point will receive a coverage value of 1,
              i.e. point has coverage.

              "opt" means optimal beam crossing angle.
              "sub" means suboptimal beam crossing angle.
              "bg" means background coverage.
              "fil" means filter coverage.
              cvg_bg, cvg_fil, and sseq_trip do not require a radar count.
              (Beam crossing angle is meaningless with one radar,
               so there is no opt or sub)
              If this file is being used, coverage values must be provided
              for all combinations of radars. """
              'cvg_opt_bg': [1, 1, 0],
              'cvg_sub_bg': [1, 1, 0],
              'cvg_opt_fil': [0, 0, 0],
              'cvg_sub_fil': [1, 1, 0],
              'cvg_bg': [1, 1, 0],
              'cvg_fil': [0, 0, 0],
              'sseq_trip': [1.0, 1.0, 0.0]}
        dda_file_name = (time_procedures.out_data_path +
                         '/dda_files/cpol_test' + year_str + month_str +
                         day_str + hour_str + minute_str + '.dda')

        pf = multidop.parameters.ParamFile(pd, dda_file_name)
        pf = multidop.parameters.CalcParamFile(pd, calc_file_name)

        # Unfortunately, text output from the analysis engine (DDA)
        # will not display until after the program completes.
        # Expect this step to take several minutes.
        bt = time.time()
        multidop.execute.run_command('./DDA ' + dda_file_name)
        print((time.time()-bt)/60.0, 'minutes to process')

        # Baseline output is not CF or Py-ART compliant.
        # This function fixes that.
        # This is why we wrote the original output to a tempfile
        # that can be safely removed.
        # The final grid will have all wind solutions outside
        # the coverage region masked.
        try:
            final_grid = multidop.grid_io.make_new_grid([grid_cpol, grid_Berr],
                                                        localfile.name)
            final_grid.write(fname)
            localfile.close()
        except:
            print('Failed to write final grid!')
            return
    else:
        print('DDA grid already exists...skipping.')
