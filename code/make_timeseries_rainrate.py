from matplotlib import use
use('agg')
import pyart
import numpy as np
import glob
import dask.bag as db
import xarray
import time
import sys
import warnings
import pandas
import os
import gc
import shutil

from distributed import Client, LocalCluster
from dask import delayed
from datetime import datetime
from netCDF4 import Dataset

warnings.filterwarnings("ignore")

excluded_fields = ['temperature', 'specific attenuation_differential_reflectivity',
                   'D0', 'NW', 'velocity_texture', 'total_power', 
                   'bringi_specific_differential_phase', 
                   'bringi_differential_phase', 'signal_to_noise_ratio', 'spectrum_width']
def get_rain_rate_over_spots(file, dis_latlon, 
                             profiler_latlon):
    
    
    rain_dis = np.nan*np.ma.ones(17)
    rain_prof = np.nan*np.ma.ones(17)
    kdp_dis = np.nan*np.ma.ones(17)
    echo_dis = np.nan*np.ma.ones(17)
    Ar_dis = np.nan*np.ma.ones(17)
    Zdr_dis = np.nan*np.ma.ones(17)
    vel_dis = np.nan*np.ma.ones(17)
    rhohv_dis = np.nan*np.ma.ones(17)
    gate_lat_dis = np.nan*np.ma.ones(17)
    gate_lon_dis = np.nan*np.ma.ones(17)
    gate_alt_dis = np.nan*np.ma.ones(17)
    kdp_prof = np.nan*np.ma.ones(17)
    echo_prof = np.nan*np.ma.ones(17)
    Ar_prof = np.nan*np.ma.ones(17)
    Zdr_prof = np.nan*np.ma.ones(17)
    rhohv_prof = np.nan*np.ma.ones(17)
    vel_prof = np.nan*np.ma.ones(17)
    gate_lat_prof = np.nan*np.ma.ones(17)
    gate_lon_prof = np.nan*np.ma.ones(17)
    gate_alt_prof = np.nan*np.ma.ones(17)
    Z_dis = np.nan*np.ma.ones(17)
    Z_prof = np.nan*np.ma.ones(17)
    Zdr_ave = np.nan*np.ma.ones(17)
    try:
        radar = pyart.io.read(file, exclude_fields=excluded_fields, 
            delay_field_loading=True)
        #atten_meta = radar.fields['specific_attenuation_reflectivity']['data']
        # Do adjusted specific attenuation retrieval
        rr = radar.fields['radar_estimated_rain_rate']
        echo = radar.fields['radar_echo_classification']
        kdp = radar.fields['corrected_specific_differential_phase']
        Ar = radar.fields['specific_attenuation_reflectivity']
        Zdr = radar.fields['corrected_differential_reflectivity']
        Z = radar.fields['reflectivity']
        vel = radar.fields['velocity']
        gate_lat = radar.gate_latitude['data']
        gate_lon = radar.gate_longitude['data']
        gate_alt = radar.gate_altitude['data']
        if('cross_correlation_ratio' in radar.fields.keys()):
            rhohv = radar.fields['cross_correlation_ratio']
        
        for i in range(radar.nsweeps):
            slice_start, slice_end = radar.get_start_end(i)
            gate_longitude = gate_lon[slice_start:slice_end,:]
            gate_latitude = gate_lat[slice_start:slice_end,:]
            gate_altitude = gate_alt[slice_start:slice_end,:]
            dist = np.sqrt((gate_latitude-dis_latlon[0])**2 + 
                           (gate_longitude-dis_latlon[1])**2)
            disdrometer_index = np.where(dist == np.min(dist))
            dist = np.sqrt((gate_latitude-profiler_latlon[0])**2 +
                           (gate_longitude-profiler_latlon[1])**2)      
            profiler_index = np.where(dist == np.min(dist))
            
            rain_dis[i] = rr['data'][slice_start+disdrometer_index[0][0], disdrometer_index[1][0]]
            rain_prof[i] = rr['data'][slice_start+profiler_index[0][0], profiler_index[1][0]]

            echo_prof[i] = echo['data'][slice_start+profiler_index[0][0], profiler_index[1][0]]
            kdp_prof[i] = kdp['data'][slice_start+profiler_index[0][0], profiler_index[1][0]]
            Ar_prof[i] = Ar['data'][slice_start+profiler_index[0][0], profiler_index[1][0]]
            Zdr_prof[i] = Zdr['data'][slice_start+profiler_index[0][0], profiler_index[1][0]]
            Z_prof[i] = Z['data'][slice_start+profiler_index[0][0], profiler_index[1][0]]
            vel_prof[i] = vel['data'][slice_start+profiler_index[0][0], profiler_index[1][0]]
            if('cross_correlation_ratio' in radar.fields.keys()):
                rhohv_prof[i] = rhohv['data'][slice_start+profiler_index[0][0], profiler_index[1][0]]

            gate_lat_prof[i] = gate_latitude[profiler_index[0][0], 
                profiler_index[1][0]]
            gate_lon_prof[i] = gate_longitude[profiler_index[0][0],
                profiler_index[1][0]]
            gate_alt_prof[i] = gate_altitude[profiler_index[0][0],
                profiler_index[1][0]]
        
            echo_dis[i] = echo['data'][slice_start+disdrometer_index[0][0], disdrometer_index[1][0]]
            kdp_dis[i] = kdp['data'][slice_start+disdrometer_index[0][0], disdrometer_index[1][0]]
            Ar_dis[i] = Ar['data'][slice_start+disdrometer_index[0][0], disdrometer_index[1][0]]
            Zdr_dis[i] = Zdr['data'][slice_start+disdrometer_index[0][0], disdrometer_index[1][0]]
             
            Zdr_vals = np.ma.array([Zdr['data'][slice_start+disdrometer_index[0][0], 
                                                disdrometer_index[1][0]],
                                    Zdr['data'][slice_start+disdrometer_index[0][0]+1,
                                                disdrometer_index[1][0]],
                                    Zdr['data'][slice_start+disdrometer_index[0][0]-1,
                                                disdrometer_index[1][0]],
                                    Zdr['data'][slice_start+disdrometer_index[0][0],
                                                disdrometer_index[1][0]+1],
                                    Zdr['data'][slice_start+disdrometer_index[0][0],
                                                disdrometer_index[1][0]-1]])
            Zdr_ave[i] = Zdr_vals.mean()
            Z_dis[i] = Z['data'][slice_start+disdrometer_index[0][0], disdrometer_index[1][0]]
            vel_dis[i] = vel['data'][slice_start+disdrometer_index[0][0], disdrometer_index[1][0]]
            if('cross_correlation_ratio' in radar.fields.keys()):
                rhohv_dis[i] = rhohv['data'][slice_start+disdrometer_index[0][0], disdrometer_index[1][0]]
 
            gate_lat_dis[i] = gate_latitude[disdrometer_index[0][0],         
                disdrometer_index[1][0]]
            gate_lon_dis[i] = gate_longitude[disdrometer_index[0][0],
                disdrometer_index[1][0]]
            gate_alt_dis[i] = gate_altitude[disdrometer_index[0][0],
                disdrometer_index[1][0]]
        
        
        print(file + ' successful!')
        
        del radar, gate_longitude, gate_latitude, gate_altitude, Z, Ar, kdp, echo, vel, Zdr, rr
        del gate_lon, gate_lat, gate_alt
        gc.collect()
    except (IOError, ValueError, KeyError, TypeError) as e:
        print(file + 'is corrupt...skipping!')
        
                       
    return np.stack([rain_dis, echo_dis, kdp_dis, Ar_dis, Zdr_dis, Z_dis,
                     vel_dis, rhohv_dis, gate_lat_dis, gate_lon_dis, gate_alt_dis, 
                     rain_prof, echo_prof, kdp_prof, Ar_prof, Zdr_prof, Z_prof,
                     vel_prof, rhohv_prof, gate_lat_prof, gate_lon_prof, gate_alt_prof,
                     Zdr_ave])


def parse_time(file):
    return datetime.strptime(file[-62:-49], '%Y%m%d_%H%M')


def ten_minute_average_impact(inp_time):
    disdrometer_path = '/lcrc/group/earthscience/rjackson/disdrometer/'
    file_name = ('twpdisdrometerC3.b1.' + "%04d" % inp_time.year +
                 "%02d" % inp_time.month + "%02d" % inp_time.day +
                 ".*.cdf")
    print('## Impact disdrometer ' + str(inp_time))
    try:
        ds = xarray.open_mfdataset(disdrometer_path+file_name)
        rain_value = ds['rain_rate'].sel(time=inp_time).values
        ds.close()
        del ds 
    except (IndexError, IOError, KeyError, ValueError, RuntimeError) as e:
        print(str(inp_time) + ' failed!')
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        return np.array([np.nan])
    return np.array([rain_value])  


def ten_minute_average_vdis(inp_time):
    disdrometer_path = '/lcrc/group/earthscience/rjackson/disdrometer/'
    file_name = ('twpvdisC3.b1.' + "%04d" % inp_time.year +
                 "%02d" % inp_time.month + "%02d" % inp_time.day +
                 ".*.cdf")
    print('## VDIS: ' + str(inp_time))
    try:
        ds = xarray.open_mfdataset(disdrometer_path+file_name)
        rain_value = ds['rain_rate'].sel(time=inp_time).values
        D0 = ds['median_volume_diameter'].sel(time=inp_time).values
        print(str(inp_time) + ' sucessfully loaded!')
        ds.close()
        del ds
    except (IOError, KeyError, ValueError, RuntimeError):
        return np.array([np.nan, np.nan])
    return np.array([rain_value, D0])


def ten_minute_average_raingauge(inp_time):
    disdrometer_path = '/lcrc/group/earthscience/rjackson/disdrometer/'
    file_name = ('twpvdisC3.b1.' + "%04d" % inp_time.year +
                 "%02d" % inp_time.month + "%02d" % inp_time.day +
                 ".*.cdf")
    print('## Rain gauge: ' + str(inp_time))
    try:
        ds = xarray.open_mfdataset(disdrometer_path+file_name)
        rain_value = ds['precip_rate'].sel(time=inp_time).values
        print(str(inp_time) + ' sucessfully loaded!')
        ds.close()
        del ds
    except (IOError, KeyError, ValueError):
        return np.nan
    return rain_value

    
def get_dros_class(inp_time):
    file_path = '/home/rjackson/data/Drosdowsky.cdf'
    in_netcdf = Dataset(file_path)
    year = in_netcdf.variables['year'][:]
    month = in_netcdf.variables['month'][:]
    day = in_netcdf.variables['day'][:]
    groups = in_netcdf.variables['groups'][:]
    the_index = np.where(np.logical_and.reduce((
        year == inp_time.year, month == inp_time.month, day == inp_time.day)))
    in_netcdf.close()
    if(the_index[0]):
        return groups[the_index[0][0]]
    else:
        return np.nan

    
def get_mjo_index(inp_time):
    mjo_index_file = '/home/rjackson/data/rmm.74toRealtime.txt'
    data = pandas.read_csv(mjo_index_file,
                           header=2,
                           delim_whitespace=True)
    data_matrix = np.ma.array(data.values)
    yearm = data_matrix[:, 0]
    monthm = data_matrix[:, 1]
    daym = data_matrix[:, 2]
    index = data_matrix[:, 5]
    the_index = np.where(np.logical_and.reduce((
        yearm == inp_time.year, monthm == inp_time.month, daym == inp_time.day)))
    if(the_index[0]):
        return index[the_index[0][0]]
    else:
        return np.nan


# Get information about disdrometer and profiler indicies in advance to speed
# up algorithm
if __name__ =='__main__':
    year = sys.argv[1]
    month = sys.argv[2]
    file_path = ('/lcrc/group/earthscience/radar/CPOL_level_1b/' +
                 '/PPI/' + year + '/' + year + month + '*' )
    output_file_path = '/lcrc/group/earthscience/rjackson/rainfall/'
    scheduler_file = '/home/rjackson/scheduler.json'
    dis_lat = -12.0 - 25/60.0 - 30/3600.0
    dis_lon = 130.0 + 53.0/60.0 + 31.2/3600.0
    profiler_lat = -12.44
    profiler_lon = 130.96
    #client = Client(scheduler_file=scheduler_file)
    cluster_path = ('cluster_' + year + month)
    if not os.path.isdir(cluster_path):
        os.mkdir(cluster_path)
    cluster = LocalCluster(n_workers=8, local_dir=cluster_path)
    client = Client(cluster)
    file_list = glob.glob(file_path + '/**/*.nc', recursive=True)
    dis_ind = (dis_lat,dis_lon)
    prof_ind = (profiler_lat, profiler_lon)
    
    bt = time.time()
    print('## Parsing dates and times...')
    file_list = sorted(file_list)
    file_bag = db.from_sequence(file_list)
    time_list = file_bag.map(parse_time).compute()

    print('## Getting Drosdowsky classification...')
    time_bag = db.from_sequence(time_list)
    dros_class = time_bag.map(get_dros_class).compute()
    dros_class = np.array(dros_class)
    
    mjo_index = time_bag.map(get_mjo_index).compute()
    mjo_index = np.array(mjo_index)

    print('## First, aggregating the video disdrometer observations to 10 min')   
    vdis_rainfall = time_bag.map(ten_minute_average_vdis).compute()
    vdis_rainfall = np.stack(vdis_rainfall)
    vdis_D0 = vdis_rainfall[:,1]
    vdis_rainfall = vdis_rainfall[:,0]

    print('## Aggregating impact disdrometer data')
    ten_minute_average_impact(time_list[0])
    dis_rainfall = time_bag.map(ten_minute_average_impact).compute()
    dis_rainfall = np.array(dis_rainfall)
    dis_rainfall = dis_rainfall[:,0]

    print('## Aggregating rain gauge data...')
    gauge_rainfall = time_bag.map(ten_minute_average_raingauge).compute()
    gauge_rainfall = np.array(gauge_rainfall)
    
    print('## Getting rainfall values...')
    dis_values = []
    #for i in range(0,len(file_list), 500):
    file_bag = db.from_sequence(file_list) 
    #the_list = file_bag.map(lambda x: get_rain_rate_over_spots(
    #                       x, dis_ind, prof_ind)).compute()
    #the_list = [get_rain_rate_over_spots(x, dis_ind, prof_ind) for x in file_list]
    #dis_values.append(np.stack(the_list, axis=0))
    #    gc.collect()
    the_futures = client.map(lambda x: get_rain_rate_over_spots(
                                  x, dis_ind, prof_ind), file_list)
    dis_values = client.gather(the_futures)
  
    for x in dis_values:
        print(x.shape)
    dis_values = np.stack(dis_values, axis=0)

    print(dis_values.shape)
    dis_rain = dis_values[:,0,:]
    echo_dis = dis_values[:,1,:]
    kdp_dis = dis_values[:,2,:]
    Ar_dis = dis_values[:,3,:]
    Zdr_dis = dis_values[:,4,:]
    Z_dis = dis_values[:,5,:]
    vel_dis = dis_values[:,6,:]
    rhohv_dis = dis_values[:,7,:]
    gate_lat_dis = dis_values[:,8,:]
    gate_lon_dis = dis_values[:,9,:]
    gate_alt_dis = dis_values[:,10,:]
    
    prof_rain = dis_values[:,11,:]
    echo_prof = dis_values[:,12,:]
    kdp_prof = dis_values[:,13,:]
    Ar_prof = dis_values[:,14,:]
    Zdr_prof = dis_values[:,15,:]
    Z_prof = dis_values[:,16,:]
    vel_prof = dis_values[:,17,:]
    rhohv_prof = dis_values[:,18,:]
    gate_lat_prof = dis_values[:,19,:]
    gate_lon_prof = dis_values[:,20,:]
    gate_alt_prof = dis_values[:,21,:]
    Zdr_ave = dis_values[:,22,:]
    print(len(time_list))
    print(dis_rain.shape)
    
    ds = xarray.Dataset({'rainfall_cpol_over_impact': (['time', 'sweep'], dis_rain),
                         'rainfall_cpol_over_vdis': (['time', 'sweep'], prof_rain),
                         'rainfall_rain_gauge': ('time', gauge_rainfall),
                         'rainfall_vdis': ('time', vdis_rainfall),
                         'rainfall_disdro': ('time', dis_rainfall),
                         'echo_dis': (['time', 'sweep'], echo_dis),
                         'echo_prof': (['time', 'sweep'], echo_prof),
                         'dros_class': (['time'], dros_class),
                         'mjo_index': (['time'], mjo_index),
                         'kdp_dis': (['time', 'sweep'], kdp_dis),
                         'Ar_dis': (['time', 'sweep'], Ar_dis),
                         'Zdr_dis': (['time', 'sweep'], Zdr_dis),
                         'Z_dis': (['time', 'sweep'], Z_dis),
                         'rhohv_dis': (['time', 'sweep'], rhohv_dis),
                         'vel_dis': (['time', 'sweep'], vel_dis),
                         'gate_lat_dis': (['time', 'sweep'], gate_lat_dis),
                         'gate_lon_dis': (['time', 'sweep'], gate_lon_dis),
                         'gate_alt_dis': (['time', 'sweep'], gate_alt_dis),
                         'kdp_prof': (['time', 'sweep'], kdp_prof),
                         'Ar_prof': (['time', 'sweep'], Ar_prof),
                         'Zdr_prof': (['time', 'sweep'], Zdr_prof),
                         'Z_prof': (['time', 'sweep'], Z_prof),
                         'rhohv_prof': (['time', 'sweep'], rhohv_prof),
                         'vel_prof': (['time', 'sweep'], vel_prof),
                         'gate_lat_prof': (['time', 'sweep'], gate_lat_prof),
                         'gate_lon_prof': (['time', 'sweep'], gate_lon_prof),
                         'gate_alt_prof': (['time', 'sweep'], gate_alt_prof),
                         'Zdr_ave': (['time', 'sweep'], Zdr_ave) 
                         },
                         coords={'time': time_list, 'sweep': np.arange(0,17)},
                         attrs={'units': 'mm hr-1', 'long_name':
                                ('Rainfall rate algorithm based on Thompson' +
                                 'et al. 2016.')})
    print(ds)
    ds.to_netcdf(path=(output_file_path + 'rainfall_rate_timeseries' + year + month + '.cdf'), mode='w')
    shutil.rmtree(cluster_path)
    client.shutdown()
