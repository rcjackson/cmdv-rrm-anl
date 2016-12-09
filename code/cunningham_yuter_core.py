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
import math
from ipyparallel import Client
from time import sleep
import time_procedures

def cunningham_yuter_grids(frame_time):
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
    from netCDF4 import Dataset
    from copy import deepcopy
    from IPython.display import Image, display
    import math
    from datetime import timedelta, datetime
    import time_procedures
    import os
    import numpy as np

    year_str = "%04d" % frame_time.year
    day_str = "%02d" % frame_time.day
    month_str = "%02d" % frame_time.month 
    hour_str = "%02d" % frame_time.hour
    minute_str = "%02d" % frame_time.minute  
    fname = ('/home/rjackson/conv_stratiform/cpol_conv_strat' + year_str 
                                                              + month_str
                                                              + day_str
                                                              + hour_str
                                                              + minute_str 
                                                              + '.nc')
    print('Does file ' + fname + ' exist?')
    if(not os.path.isfile(fname)):
        Radar = time_procedures.get_radar_from_cpol(frame_time)

        # Calculate texture of radial velocity field and mask reflectivity field
        print('Calculating texture....')
        nyq_Gunn = Radar.instrument_parameters['nyquist_velocity']['data'][0]

        data = ndimage.filters.generic_filter(Radar.fields['velocity']['data'],
                                              pyart.util.interval_std, size = (6,6),
                                              extra_arguments = (-nyq_Gunn, nyq_Gunn))
        filtered_data = ndimage.filters.median_filter(data, size = (4,4))
        texture_field = pyart.config.get_metadata('velocity')
        texture_field['data'] = filtered_data
        cp = deepcopy(Radar.fields['corrected_reflectivity']['data'])
        texture = texture_field['data']
        Radar.add_field('velocity_texture', texture_field, replace_existing = True)    
	cp = np.ma.masked_where(texture > 4, cp)
        Radar.add_field_like('corrected_reflectivity', 
                             'reflectivity', 
                             cp,
                             replace_existing=True)
	
	
        def grid_radar(radar, grid_shape=(20, 301, 301), xlim=(-150000, 150000),
               ylim=(-150000, 150000), zlim=(1000, 20000), bsp=1.0, 
               min_radius=750, h_factor=4.0, nb=1.5,
               fields=['DT', 'VT'], origin=None):

            bt = time.time()
            radar_list = [radar]
            if origin is None:
                origin = (radar.latitude['data'][0],
                          radar.longitude['data'][0])
            grid = pyart.map.grid_from_radars(
                                              radar_list, grid_shape=grid_shape,
                                              grid_limits=(zlim, ylim, xlim),
                                              grid_origin=origin, fields=fields,
                                              weighting_function='Cressman',
                                              gridding_algo='map_gates_to_grid',
                                              h_factor=h_factor,
                                              min_radius=min_radius,
                                              bsp=bsp,
                                              nb=nb)
            print(time.time() - bt, 'seconds to grid radar')
            return grid
	
        grid_cpol = grid_radar(Radar, 
                               origin=(Radar.latitude['data'][0], Radar.longitude['data'][0]),
                               xlim=(-150000, 150000), ylim=(-150000, 150000), 
                               fields=['reflectivity'], min_radius=750.0, bsp=1.0, nb=1.5,
                               h_factor=2.0,
                               zlim=(500, 20000), grid_shape=(40, 121, 121))
	
        convective = pyart.retrieve.steiner_conv_strat(grid_cpol) 

        # Write convective-stratiform classification into netCDF file
        out_netcdf = Dataset(fname, mode='w')            
        out_netcdf.createDimension('X', len(grid_cpol.x['data']))
        out_netcdf.createDimension('Y', len(grid_cpol.y['data']))

        x_file = out_netcdf.createVariable('Grid_X', 'f4', ('X',))
        x_file.long_name = 'Distance north of CPOL'
        x_file.units = 'km'
        x_file[:] = grid_cpol.x['data']/1e3
        
        y_file = out_netcdf.createVariable('Grid_Y', 'f4', ('Y',))
        y_file.long_name = 'Distance east of CPOL'
        y_file.units = 'km'
        y_file[:] = grid_cpol.y['data']/1e3
         
        conv_file = out_netcdf.createVariable('strat_conv', 'i4', ('X','Y'))
        conv_file.long_name = 'Steiner convective classification'
        conv_file.units = '0 = Little precip, 1 = Stratiform, 2 = Convective'
        conv_file[:] = convective['data']

        out_netcdf.close()
    else:
        print('Classification already exists...skipping.')


