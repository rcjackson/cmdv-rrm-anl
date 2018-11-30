from matplotlib import use
use('agg')
import pydsd
import dask.bag as db
import xarray
import numpy as np
import pytmatrix
import datetime
import pandas

from glob import glob
from datetime import datetime
from netCDF4 import Dataset

jw_disdrometer_file_path = '/lcrc/group/earthscience/rjackson/disdrometer/twpvdis*.cdf'
dros_file_path = '/home/rjackson/data/Drosdowsky.cdf'
mjo_index_file = '/home/rjackson/data/rmm.74toRealtime.txt'

dros_cdf = Dataset(dros_file_path)
year = dros_cdf.variables['year'][:]
month = dros_cdf.variables['month'][:]
day = dros_cdf.variables['day'][:]
groups = dros_cdf.variables['groups'][:]

data = pandas.read_csv(mjo_index_file,
                           header=2,
                           delim_whitespace=True)
data_matrix = np.ma.array(data.values)
yearm = data_matrix[:, 0]
monthm = data_matrix[:, 1]
daym = data_matrix[:, 2]
index = data_matrix[:, 5]


def get_PSD(filename):
    def brandes2002(D_eq):
        return 0.9951 + 0.0251*D_eq - 0.03644*D_eq**2 + 0.005303*D_eq**3 - 0.0002492*D_eq**4

    #try:
    JWReader = pydsd.read_arm_vdis_b1(filename)
    JWReader.set_scattering_temperature_and_frequency(scattering_temp=20,
                                                      scattering_freq=5.61e9)
    JWReader.calculate_radar_parameters(dsr_func=brandes2002, std_c=12.0)
    #except KeyError:
    #    print('No nd key!')
    #    return []
    print(filename + ' processed!')
    return JWReader


def get_dros_class(inp_time):
    the_index = np.where(np.logical_and.reduce((
        year == inp_time.year, month == inp_time.month, day == inp_time.day)))
    if(the_index[0]):
        return groups[the_index[0][0]]
    else:
        return np.nan

    
def get_mjo_index(inp_time):
    the_index = np.where(np.logical_and.reduce((
        yearm == inp_time.year, monthm == inp_time.month, daym == inp_time.day)))
    if(the_index[0]):
        return index[the_index[0][0]]
    else:
        return np.nan


if __name__ == '__main__':
    # Get radar parameters, DSDs fron entire JW Disdrometer dataset
    file_list = sorted(glob(jw_disdrometer_file_path))    
    file_bag = db.from_sequence(file_list)
    print('Doing disdrometer calculations...')    
    the_PSDs = file_bag.map(get_PSD).compute()
    the_PSDs = [x for x in the_PSDs if x != []]
    print(the_PSDs[0].fields.keys())
    Ai = np.ma.concatenate([x.fields['Ai']['data'] for x in the_PSDs])
    rain_rate = np.ma.concatenate([x.fields['rain_rate']['data'] for x in the_PSDs])
    Nd = np.ma.concatenate([x.fields['Nd']['data'] for x in the_PSDs], axis=0)
    Zh = np.ma.concatenate([x.fields['Zh']['data'] for x in the_PSDs])
    Kdp = np.ma.concatenate([x.fields['Kdp']['data'] for x in the_PSDs]) 
    num_drop = np.ma.concatenate([x.fields['Nd']['data'] for x in the_PSDs])
    #liq_water = np.ma.concatenate([x.fields['liq_water']['data'] for x in the_PSDs])
    n_0 = np.ma.concatenate([x.fields['n_0']['data'] for x in the_PSDs])
    #d_max = np.ma.concatenate([x.fields['d_max']['data'] for x in the_PSDs])
    Adr = np.ma.concatenate([x.fields['Adr']['data'] for x in the_PSDs])
    Zdr = np.ma.concatenate([x.fields['Zdr']['data'] for x in the_PSDs])
    times = []
    times_dt = [] 
    for i in range(len(the_PSDs)):
        y = [np.datetime64(datetime.utcfromtimestamp(x)) for x in the_PSDs[i].time['data']]
        y = np.array(y, dtype='datetime64[ns]')
        dt = np.array([datetime.utcfromtimestamp(x) for x in the_PSDs[i].time['data']])
        if(len(y) > 0):
            times.append(y)
            times_dt.append(dt)
    times_dt = np.concatenate(times_dt)
    time_bag = db.from_sequence(times_dt)
    dros_class = time_bag.map(get_dros_class).compute() 
    dros_class = np.array(dros_class)
    mjo_index = time_bag.map(get_mjo_index).compute()
    mjo_index = np.array(mjo_index)            
    times = np.concatenate(times)
    bin_edges = the_PSDs[0].bin_edges['data']
    bin_edges = (bin_edges[:-1]+bin_edges[1:])/2 
    
    ds = xarray.Dataset({'Ai': (['time'], Ai),
                         'rain_rate': (['time'], rain_rate), 
                         'Nd': (['time', 'bin_edges'], Nd),
                         'Zh': (['time'], Zh),
                         'Kdp': (['time'], Kdp),
                         'num_drop': (['time', 'bin_edges'], num_drop),
                         'n_0': (['time'], n_0),
                         'Adr': (['time'], Adr),
                         'Zdr': (['time'], Zdr),
                         'dros_regime': (['time'], dros_class),
                         'mjo_index': (['time'], mjo_index)},
                        coords={'time': times, 'bin_edges': bin_edges})
    print(ds)
    ds.to_netcdf('PSDs_w_radar_params_vdis.cdf')
