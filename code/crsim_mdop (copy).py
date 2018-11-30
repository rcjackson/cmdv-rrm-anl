import multidop
import tempfile
import pyart
import sys
import numpy as np
import dask.bag as db
import os
import time

from distributed import Client, wait, LocalCluster
from sklearn.metrics import mean_squared_error
from netCDF4 import Dataset

cpol_file_path = '/home/rjackson/runtwp1a_cfradial/crsimwrfruntwp1acpolCPOL.a02c.20060120.005000.cdf'
berr_file_path = '/home/rjackson/runtwp1a_cfradial/crsimwrfruntwp1acpolBERR.a02c.20060120.005000.cdf'
grid_cpol_path = '/home/rjackson/runtwp1a_cfradial/berr_Darwin.nc'
grid_berr_path = '/home/rjackson/runtwp1a_cfradial/cpol_Darwin.nc'
dda_file_path = '/lcrc/group/earthscience/rjackson/crsim_mdop/dda/'
frprmn_file_path = '/lcrc/group/earthscience/rjackson/crsim_mdop/frprmn/'
out_grid_path = '/lcrc/group/earthscience/rjackson/crsim_mdop/grids/'
wrf_run_path = '/lcrc/group/earthscience/rjackson/satoshi_runs/runtwp1a/wrfout_d01_2006-01-20_00:00:00'


def do_test_retrieval(const_tuple, origin_lon, origin_lat, grid1, grid2):
    data_weighting_constant = const_tuple[0]
    mass_continuity_constant = const_tuple[1]
    vorticity_constant = const_tuple[2]
    smoothing_constant = const_tuple[3]

    dx = grid1.x['data'][1]-grid1.x['data'][0]
    dy = grid1.y['data'][1]-grid1.y['data'][0]
    dz = grid1.z['data'][1]-grid1.z['data'][0]
    nx = len(grid1.x['data'])
    ny = len(grid1.y['data'])
    nz = len(grid1.z['data'])
    
    # You don't have to define everything. Most of these keywords are default values.
    # If you don't define something the program will provide a default value.
    # Check parameters.py for what keyword default values are.
    localfile = tempfile.NamedTemporaryFile()
    fname = (frprmn_file_path + 'frprmn' + str(mass_continuity_constant) +
             'v' + str(vorticity_constant) + 
             's' + str(smoothing_constant) + '.nc')
    pd = {'dir': './',
          'x': [grid1.x['data'][0], dx, nx],   # start, step, max = min + (steps-1)
          'y': [grid1.y['data'][0], dy, ny],
          'z': [grid1.z['data'][0], dz, nz],
          'grid': [origin_lon, origin_lat, 50.0],
          'files': [grid_cpol_path,
                    grid_berr_path],
          'radar_names': ['Berrima', 'CPOL'],
          'refl': 'DT',  # Name of reflectivity field. Must be common between radars.
          'vt': 'VT',  # Name of velocity field. Must be common between radars.
          'bgfile': 'sounding_file', # Name of sounding file
          'writeout': localfile.name, # Name of output grid file
          'frprmn_out': frprmn_out_name,
          'min_cba': 30.0,  # Minimum beam-crossing angle
          'calc_params': 'cpol_calc.dda', # .dda file for parameters related to minimization routine
          'anel': 1, # 0 = Boussinesq approximation for mass conservation, 1 = anelastic 
          'laplace': 0, # 0 = 1st order derivatives for smoothing, 1 = second
          'read_dataweights': 2, # 0 = calculate data constraint weights/output, 1 = read from file, 2 = weigh all equally
          'max_dist': 10.0, # How much distance analysis and observational grid must match in m
          'cutoff': 0.0, # Deny observations below this level from analysis (m)
          'UT': 0.0, # U of prescribed storm motion vector
          'VT': 0.0, # V of prescribed storm motion vector
          'output_error': 0, # 1 = output verification stats after each iteration
          'weak_height': -1, # Sounding height constraint weakened in regions > 10 dBZ below this height (-1 = disabled)
          'upper_bc': 1, # 1 = w = 0 as upper boundary condition, -1 = ignore
          'itmax_frprmn': [200, 10], # max iterations in frprmn function
          'itmax_dbrent': 200, # max iterations in dbrent function
          'C1b': data_weighting_constant,  # Data weighting factor
          'C2b': mass_continuity_constant,  # Mass continuity weighting factor
          'C3b': vorticity_constant,  # Vorticity weighting factor
          'C4b': smoothing_constant,  # Horizontal smoothing factor
          'C5b': smoothing_constant,  # Vertical smoothing factor
          'C8b': 0.000,  # Sounding factor
          'vary_weights': 0,
          # Define filter with ONE of the following forms.
          # filter: None,
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
    fname = (out_grid_path + 'mdop_testm' + str(mass_continuity_constant) +
             'v' + str(vorticity_constant) + 
             's' + str(smoothing_constant) + '.nc')
    dda_file_name = (dda_file_path + 'cpol_testm' + str(mass_continuity_constant) +
        'v' + str(vorticity_constant) + 's' + str(smoothing_constant) + '.dda')
    param_file_name =  (dda_file_path + 'cpol_calcm' + str(mass_continuity_constant) +
        'v' + str(vorticity_constant) + 's' + str(smoothing_constant) + '.dda')
    if(os.path.isfile(fname)):
        print('Grid already exists! Skipping grid generation.')
    else:
        pf = multidop.parameters.ParamFile(pd, dda_file_name)
        pf = multidop.parameters.CalcParamFile(pd, param_file_name)
        # Unfortunately, text output from the analysis engine (DDA) will not display
        # until after the program completes. Expect this step to take several minutes.
        bt = time.time()
        multidop.execute.do_analysis(dda_file_name, 
                                    cmd_path='/home/rjackson/multidop/src/DDA')
        print((time.time()-bt)/60.0, 'minutes to process')
        # Baseline output is not CF or Py-ART compliant. This function fixes that.
        # This is why we wrote the original output to a tempfile that can be safely removed.
        # The final grid will have all wind solutions outside the coverage region masked.
        final_grid = multidop.grid_io.make_new_grid([grid1, grid2], 
                                                    localfile.name,
                                                    use_mask=True,
                                                    min_refl=-40.0)
        final_grid.write(fname)
        localfile.close()
      

# Default scale @ 1, 1 km horiz, 0.5 km vertical
def regrid_wrf_data(timestep, scale_horiz, scale_vert):
    # Load wrf data
    wrf_cdf = Dataset(wrf_run_path, mode='r')
    Z_wrf = wrf_cdf.variables['REFL_10CM'][:,:,:,:]
    Lat_wrf = wrf_cdf.variables['XLAT'][:,:,:]
    Lon_wrf = wrf_cdf.variables['XLONG'][:,:,:]
    W_wrf = wrf_cdf.variables['W'][:]
    V_wrf = wrf_cdf.variables['V'][:]
    U_wrf = wrf_cdf.variables['U'][:]
    PH_wrf = wrf_cdf.variables['PH'][:]
    PHB_wrf = wrf_cdf.variables['PHB'][:]
    ETA_wrf = wrf_cdf.variables['ZNW'][:]
    array_shape = PH_wrf.shape
    alt_wrf = (PH_wrf+PHB_wrf)/9.81
    
    w = W_wrf[timestep,:-1,:,:]
    v = V_wrf[timestep,:,:-1,:]
    u = U_wrf[timestep,:,:,:-1]
    z = Z_wrf[timestep,:,:,:]
    cpol_x = 56-87
    cpol_y = 100-87
    berr_x = 56-87
    berr_y = 74-87
    levels = np.arange(0,20000,500*scale_vert)
    # Set new grid to be same as multidop grid spec
    new_grid_x, new_grid_y, new_grid_z = np.meshgrid(
            np.arange(-55.,50.,scale_horiz)+cpol_x, 
            np.arange(-50.,30.,scale_horiz)+cpol_y, 
            levels)
    x = range(-87,87,1)
    y = range(-87,87,1)
    z = alt_wrf[timestep,:-1,1,1]
    print(v.shape)
    w_interp = RegularGridInterpolator((z, y, x), w)
    v_interp = RegularGridInterpolator((z, y, x), v)
    u_interp = RegularGridInterpolator((z, y, x), u)
    z_interp = RegularGridInterpolator((z, y, x), Z_wrf[timestep,:,:,:])
    w_new = w_interp((new_grid_z, new_grid_y, new_grid_x))
    v_new = v_interp((new_grid_z, new_grid_y, new_grid_x))
    u_new = u_interp((new_grid_z, new_grid_y, new_grid_x))
    z_new = z_interp((new_grid_z, new_grid_y, new_grid_x))
    
    # Mask out entries where we are outside of optimal beam crossing angle
    bca = get_bca_wrf(new_grid_x, new_grid_y)
    for i in range(len(new_grid_z)):
        w_new[i] = np.ma.masked_where(np.logical_or(
                bca < np.pi/6, bca > 5*np.pi/6), w_new[i])
        v_new[i] = np.ma.masked_where(np.logical_or(
                bca < np.pi/6, bca > 5*np.pi/6), v_new[i])
        u_new[i] = np.ma.masked_where(np.logical_or(
                bca < np.pi/6, bca > 5*np.pi/6), u_new[i])
        z_new[i] = np.ma.masked_where(np.logical_or(
                bca < np.pi/6, bca > 5*np.pi/6), z_new[i])
    
    wrf_cdf.close()
    
    return u_new, v_new, w_new, z_new, levels 
    

def get_bca(grid):
    berr_origin = [-0.727681,-26557.7]
    x,y = np.meshgrid(grid.x['data'], grid.y['data'])
    a = np.sqrt(np.multiply(x,x)+np.multiply(y,y))
    b = np.sqrt(pow(x-berr_origin[0],2)+pow(y-berr_origin[1],2))
    c = np.sqrt(berr_origin[0]*berr_origin[0]+berr_origin[1]*berr_origin[1])
    theta_1 = np.arccos(x/a)
    theta_2 = np.arccos((x-berr_origin[1])/b)
    return np.arccos((a*a+b*b-c*c)/(2*a*b))


def get_bca_wrf(x,y):
    cpol_origin = np.array([56, 100])-87
    berr_origin = np.array([56, 74])-87-cpol_origin
    x = x-cpol_origin[0]
    y = y-cpol_origin[1]
    a = np.sqrt(np.multiply(x,x)+np.multiply(y,y))
    b = np.sqrt(pow(x-berr_origin[0],2)+pow(y-berr_origin[1],2))
    c = np.sqrt(berr_origin[0]*berr_origin[0]+berr_origin[1]*berr_origin[1])
    theta_1 = np.arccos(x/a)
    theta_2 = np.arccos((x-berr_origin[1])/b)
    return np.arccos((a*a+b*b-c*c)/(2*a*b))


def calc_rmse(const_tuple, u_new, v_new, w_new, levels):
    # Then, restrict further to DD lobes
    rmse_u = np.zeros(len(levels))
    rmse_v = np.zeros(len(levels))
    rmse_w = np.zeros(len(levels))
    
    # Load multidop grid
    data_weighting_constant = const_tuple[0]
    mass_continuity_constant = const_tuple[1]
    vorticity_constant = const_tuple[2]
    smoothing_constant = const_tuple[3]
    fname = (out_grid_path + 'mdop_testm' + str(mass_continuity_constant) +
             'v' + str(vorticity_constant) + 
             's' + str(smoothing_constant) + '.nc')
    
    the_grid = pyart.io.read_grid(fname)
    u_mdop = the_grid.fields['eastward_wind']['data']
    v_mdop = the_grid.fields['northward_wind']['data']
    w_mdop = the_grid.fields['upward_air_velocity']['data']
    
    
    bca_mdop = get_bca(the_grid)
    for i in range(len(levels)):
        # Get beam crossing angle for mask
        u_mdop[i] = np.ma.masked_where(np.logical_or(
                bca < np.pi/6, bca > 5*np.pi/6), u_mdop[i])
        v_mdop[i] = np.ma.masked_where(np.logical_or(
                bca < np.pi/6, bca > 5*np.pi/6), v_mdop[i])
        w_mdop[i] = np.ma.masked_where(np.logical_or(
                bca < np.pi/6, bca > 5*np.pi/6), w_mdop[i])
        
        # Get RMSE
        comp_indicies = np.where(np.logical_and(u_mdop.mask[i] == False,
                                                u_new.mask[i] == False))
        rmse_u[i] = np.sqrt(mean_squared_error(u_new[i], u_mdop[i]))
        comp_indicies = np.where(np.logical_and(v_mdop.mask[i] == False,
                                                v_new.mask[i] == False))
        rmse_v[i] = np.sqrt(mean_squared_error(v_new[i], v_mdop[i]))
        comp_indicies = np.where(np.logical_and(w_mdop.mask[i] == False,
                                                w_new.mask[i] == False))
        rmse_w[i] = np.sqrt(mean_squared_error(w_new[i], w_mdop[i]))
    
        
    # Close multidop grid to save memory
    del the_grid
    return np.stack([rmse_u, rmse_v, rmse_w])
