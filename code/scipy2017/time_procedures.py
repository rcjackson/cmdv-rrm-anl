# This module handles all of the time lookups for soundings and radar
# All of the file_name_str entries will have to be adjusted to fit to
# your radar dataset's naming convention

import glob
import numpy as np
import math
import matplotlib
matplotlib.use('agg')
import pyart
import time
from copy import deepcopy

data_path_berr = ('/lcrc/group/earthscience/radar/stage/radar_disk_two' +
                  '/berr_rapic/')
out_data_path = '/lcrc/group/earthscience/rjackson/multidop_grids/'
cpol_grid_data_path = '/lcrc/group/earthscience/rjackson/data/radar/grids/'
data_path_sounding = '/lcrc/group/earthscience/rjackson/soundings/'
berr_data_file_path = ('/lcrc/group/earthscience/radar/stage/' +
                       '/radar_disk_two/berr_rapic/')
data_path_cpol = ('/lcrc/group/earthscience/radar/stage/radar_disk_two/' +
                  '/cpol_rapic/')
data_path_cpol_cfradial = '/lcrc/group/earthscience/rjackson/cpol'
data_path_berr_cfradial = '/lcrc/group/earthscience/rjackson/berr'
out_file_path = '/lcrc/group/earthscience/rjackson/quicklook_plots/'


# Get a Radar object given a time period in the CPOL dataset
def get_radar_from_berr(time):
    from datetime import timedelta, datetime
    year_str = "%04d" % time.year
    month_str = "%02d" % time.month
    day_str = "%02d" % time.day
    hour_str = "%02d" % time.hour
    minute_str = "%02d" % time.minute
    second_str = "%02d" % time.second
    file_name_str = (data_path_berr +
                     'BerrimaVol' +
                     year_str +
                     month_str +
                     day_str +
                     '_' +
                     hour_str +
                     minute_str +
                     second_str +
                     '_deal.uf')
    radar = pyart.io.read(file_name_str)
    return radar


""" get_grid_times_cpol
     start_year = Start year of animation
     start_month = Start month of animation
     start_day = Start day of animation
     start_hour = Start hour of animation
     end_year = End year of animation
     end_month = End month of animation
     end_day = End day of animation
     end_minute = End minute of animation
     minute_interval = Interval in minutes between scans (default is 5)
     This procedure acquires an array of Grid classes between
     start_time and end_time  """


def get_grid_times_cpol(start_year, start_month, start_day,
                        start_hour, start_minute, end_year,
                        end_month, end_day, end_hour,
                        end_minute, minute_interval=5):

    from datetime import timedelta, datetime
    start_time = datetime(start_year,
                          start_month,
                          start_day,
                          start_hour,
                          start_minute,
                          )
    end_time = datetime(end_year,
                        end_month,
                        end_day,
                        end_hour,
                        end_minute,
                        )

    deltatime = end_time - start_time

    if(deltatime.seconds > 0 or deltatime.minute > 0):
        no_days = deltatime.days + 1
    else:
        no_days = deltatime.days

    if(start_day != end_day):
        no_days = no_days + 1

    days = range(0, no_days)
    print('We are about to load grid files for ' + str(no_days) + ' days')

    # Find the list of files for each day
    cur_time = start_time

    file_list = []
    time_list = []
    date_list_final = []
    for i in days:
        year_str = "%04d" % cur_time.year
        day_str = "%02d" % cur_time.day
        month_str = "%02d" % cur_time.month
        dir_str = year_str + '/' + month_str + '/' + day_str + '/'
        format_str = (cpol_grid_data_path + dir_str + 'cpol_' + year_str +
                      month_str + day_str + '*.nc')

        print('Looking for files with format ' + format_str)

        data_list = glob.glob(format_str)
        if(len(data_list) > 0):
            day = datetime(cur_time.year, cur_time.month, cur_time.day,
                           0, 0, 1)
            date_list_final.append(day)

        for j in range(0, len(data_list)):
            file_list.append(data_list[j])
        cur_time = cur_time + timedelta(days=1)

    # Parse all of the dates and time in the interval
    # and add them to the time list
    past_time = []
    for file_name in file_list:
        date_str = file_name[-15:-3]
        year_str = date_str[0:4]
        month_str = date_str[4:6]
        day_str = date_str[6:8]
        hour_str = date_str[8:10]
        minute_str = date_str[10:12]
        cur_time = datetime(int(year_str), int(month_str), int(day_str),
                            int(hour_str), int(minute_str), 0)
        time_list.append(cur_time)

    # Sort time list and make sure time are at least xx min apart
    time_list.sort()
    time_list_sorted = deepcopy(time_list)

    time_list_final = []
    past_time = []

    for times in time_list_sorted:
        cur_time = times

        if(past_time == []):
            past_time = cur_time

        if(cur_time - past_time >= timedelta(minutes=minute_interval) and
           cur_time >= start_time and cur_time <= end_time):
            time_list_final.append(cur_time)
            past_time = cur_time

    return time_list_final

""" get_radar_times_cpol
     start_year = Start year of animation
     start_month = Start month of animation
     start_day = Start day of animation
     start_hour = Start hour of animation
     end_year = End year of animation
     end_month = End month of animation
     end_day = End day of animation
     end_minute = End minute of animation
     minute_interval = Interval in minutes between scans (default is 5)
     This procedure acquires an array of Radar
     classes between start_time and end_time.  """


def get_radar_times_cpol(start_year, start_month, start_day,
                         start_hour, start_minute, end_year,
                         end_month, end_day, end_hour,
                         end_minute, minute_interval=1):

    from datetime import timedelta, datetime
    from parse import parse
    start_time = datetime(start_year, start_month, start_day,
                          start_hour, start_minute,)
    end_time = datetime(end_year, end_month, end_day,
                        end_hour, end_minute,)

    deltatime = end_time - start_time

    if(deltatime.seconds > 0 or deltatime.minute > 0):
        no_days = deltatime.days + 1
    else:
        no_days = deltatime.days

    if(start_day != end_day):
        no_days = no_days + 1

    days = range(0, no_days)
    print(('We are about to load grid files for ' +
           str(no_days) + ' days'))

    # Find the list of files for each day
    cur_time = start_time

    file_list = []
    time_list = []
    date_list_final = []
    for i in days:
        year_str = "%04d" % cur_time.year
        day_str = "%02d" % cur_time.day
        month_str = "%02d" % cur_time.month
        # Adjust to your dataset
        if(cur_time.year > 2007):
            format_str = (data_path_cpol_cfradial + '/' + year_str + '/' +
                          year_str + month_str + day_str + '/cfrad.' +
                          year_str + month_str + day_str + '*UNKNOWN_SUR.nc')
        else:
            format_str = (data_path_cpol_cfradial + '/' + year_str + '/' +
                          year_str + month_str + day_str + '/Gunn_pt*' +
                          year_str + month_str + day_str + '*ppi.nc')

        print('Looking for files with format ' + format_str)
        data_list = glob.glob(format_str)

        if(len(data_list) > 0):
            day = datetime(cur_time.year, cur_time.month, cur_time.day,
                           0, 0, 1)
            date_list_final.append(day)

        for j in range(0, len(data_list)):
            file_list.append(data_list[j])
        cur_time = cur_time + timedelta(days=1)

    # Parse all of the dates and time in the interval and
    # add them to the time list
    past_time = []
    for file_name in file_list:
        if(not file_name[-6:] == 'ppi.nc'):
            new_format_str = (data_path_cpol_cfradial +
                              '/{:d}/{:d}/' +
                              'cfrad.{:d}_{:d}.{:d}_to_{:d}_{:d}' +
                              '.{:d}_Gunn_Pt_v{:d}_UNKNOWN_SUR.nc')
            print(file_name)
            parameters = parse(new_format_str, file_name)
            year_str = np.floor(parameters[2]/10000)
            month_str = np.floor((parameters[2] - year_str*10000)/100)
            day_str = np.floor(parameters[2] - year_str*10000 - month_str*100)
            hour_str = np.floor(parameters[3]/10000)
            minute_str = np.floor((parameters[3] - hour_str*10000)/100)
            second_str = np.floor(parameters[3] - hour_str*10000 -
                                  minute_str*100)
        else:
            date_str = file_name[-20:-6]
            year_str = date_str[0:4]
            month_str = date_str[4:6]
            day_str = date_str[6:8]
            hour_str = date_str[8:10]
            minute_str = date_str[10:12]
            second_str = date_str[12:14]

        print(year_str)
        cur_time = datetime(int(year_str), int(month_str), int(day_str),
                            int(hour_str), int(minute_str), int(second_str))
        time_list.append(cur_time)

    # Sort time list and make sure time are at least xx min apart
    time_list.sort()
    time_list_sorted = deepcopy(time_list)

    time_list_final = []
    past_time = []

    for times in time_list_sorted:
        cur_time = times

        if(past_time == []):
            past_time = cur_time

        if(cur_time >= start_time and cur_time <= end_time):
            time_list_final.append(cur_time)
            past_time = cur_time

    return time_list_final, date_list_final

""" get_radar_times_cpol
    start_year = Start year of animation
    start_month = Start month of animation
    start_day = Start day of animation
    start_hour = Start hour of animation
    end_year = End year of animation
    end_month = End month of animation
    end_day = End day of animation
    end_minute = End minute of animation
    minute_interval = Interval in minutes between scans (default is 5)
    This procedure acquires an array of Radar classes
    between start_time and end_time  """


def get_radar_times_berr(start_year, start_month, start_day,
                         start_hour, start_minute, end_year,
                         end_month, end_day, end_hour,
                         end_minute, minute_interval=5):

    from datetime import timedelta, datetime
    from parse import parse
    start_time = datetime(start_year, start_month, start_day,
                          start_hour, start_minute,)
    end_time = datetime(end_year, end_month, end_day, end_hour, end_minute,)

    deltatime = end_time - start_time
    if(deltatime.seconds > 0 or deltatime.minute > 0):
        no_days = deltatime.days + 1
    else:
        no_days = deltatime.days

    if(start_day != end_day):
        no_days = no_days + 1

    days = range(0, no_days)
    print('We are about to load grid files for ' + str(no_days) + ' days')

    # Find the list of files for each day
    cur_time = start_time
    file_list = []
    time_list = []
    date_list_final = []
    for i in days:
        year_str = "%04d" % cur_time.year
        day_str = "%02d" % cur_time.day
        month_str = "%02d" % cur_time.month
        # Adjust to your dataset
        format_str = (data_path_berr_cfradial +
                      '/' + year_str + '/' + year_str +
                      month_str + day_str + '/cfrad.' +
                      year_str + month_str + day_str + '*.nc')

        print('Looking for files with format ' + format_str)
        data_list = glob.glob(format_str)
        if(len(data_list) > 0):
            day = datetime(cur_time.year, cur_time.month, cur_time.day,
                           0, 0, 1)
            date_list_final.append(day)

        for j in range(0, len(data_list)):
            file_list.append(data_list[j])
        cur_time = cur_time + timedelta(days=1)

    # Parse all of the dates and time in the interval
    # and add them to the time list
    past_time = []
    for file_name in file_list:
        if(file_name[-13:] == 'el0.50_SUR.nc'):
            new_format_str = (data_path_berr_cfradial +
                              '/{:d}/{:d}/' +
                              'cfrad.{:d}_{:d}.{:d}_to_{:d}_{:d}.' +
                              '{:d}_Berr_v{:d}_s{:d}_el0.50_SUR.nc')
        else:
            new_format_str = (data_path_berr_cfradial +
                              '/{:d}/{:d}/' +
                              'cfrad.{:d}_{:d}.{:d}_to_{:d}_{:d}.' +
                              '{:d}_Berrima_v{:d}_UNKNOWN_SUR.nc')

        parameters = parse(new_format_str, file_name)
        year_str = np.floor(parameters[2]/10000)
        month_str = np.floor((parameters[2] - year_str*10000)/100)
        day_str = np.floor(parameters[2] - year_str*10000 - month_str*100)
        hour_str = np.floor(parameters[3]/10000)
        minute_str = np.floor((parameters[3] - hour_str*10000)/100)
        second_str = np.floor(parameters[3] - hour_str*10000 - minute_str*100)
        cur_time = datetime(int(year_str), int(month_str), int(day_str),
                            int(hour_str), int(minute_str), int(second_str))
        # Only interested at scans at 10's of minutes since
        # they are full volumes -- adjust for your dataset
        if(int(minute_str) % 10 == 0):
            time_list.append(cur_time)

    # Sort time list and make sure time are at least xx min apart
    time_list.sort()
    time_list_sorted = deepcopy(time_list)

    time_list_final = []
    past_time = []

    for times in time_list_sorted:
        cur_time = times

        if(past_time == []):
            past_time = cur_time

        if(cur_time >= start_time and cur_time <= end_time):
            time_list_final.append(cur_time)
            past_time = cur_time

    return time_list_final, date_list_final


# Write to cfradial file given a time (useful for adding fields)
def write_radar_to_cpol(radar, time):
    import glob
    year_str = "%04d" % time.year
    month_str = "%02d" % time.month
    day_str = "%02d" % time.day
    hour_str = "%02d" % time.hour
    minute_str = "%02d" % time.minute
    second_str = "%02d" % time.second
    if(time.year > 2007):
        file_name_str = (data_path_cpol_cfradial +
                         '/' + year_str + '/' + year_str + month_str +
                         day_str + '/' + 'cfrad.' + year_str + month_str +
                         day_str + '_' + hour_str + minute_str + '*.nc')
    else:
        file_name_str = (data_path_cpol_cfradial + '/' + year_str + '/' +
                         year_str + month_str + day_str + '/' + 'Gunn_pt_' +
                         year_str + month_str + day_str + hour_str +
                         minute_str + second_str + '*.nc')

    file_name = glob.glob(file_name_str)
    file_name = file_name[0]
    pyart.io.write_cfradial(file_name, radar)


def write_radar_to_berr(radar, time):
    import glob
    year_str = "%04d" % time.year
    month_str = "%02d" % time.month
    day_str = "%02d" % time.day
    hour_str = "%02d" % time.hour
    minute_str = "%02d" % time.minute
    second_str = "%02d" % time.second
    file_name_str = (data_path_berr_cfradial +
                     '/' + year_str + '/' + year_str + month_str + day_str +
                     '/cfrad.' + year_str + month_str + day_str + '_' +
                     hour_str + minute_str + '*.nc')
    file_name = glob.glob(file_name_str)
    file_name = file_name[0]
    pyart.io.write_cfradial(file_name, radar)


# Get a Radar object given a time period in the CPOL dataset
def get_radar_from_cpol(time):
    from datetime import timedelta, datetime
    import glob
    year_str = "%04d" % time.year
    month_str = "%02d" % time.month
    day_str = "%02d" % time.day
    hour_str = "%02d" % time.hour
    minute_str = "%02d" % time.minute
    second_str = "%02d" % time.second

    if(time.year <= 2007):
        file_name_str = (data_path_cpol_cfradial +
                         '/' + year_str + '/' + year_str +
                         month_str + day_str + '/Gunn_pt_' +
                         year_str + month_str + day_str +
                         hour_str + minute_str + '*.nc')
    else:
        file_name_str = (data_path_cpol_cfradial +
                         '/' + year_str + '/' + year_str + month_str +
                         day_str + '/cfrad.' + year_str + month_str +
                         day_str + '_' + hour_str + minute_str + '*.nc')
    print('Opening ' + file_name_str)
    file_name = glob.glob(file_name_str)
    radar = pyart.io.read(file_name[0])
    return radar


# Get a Radar object given a time period in the CPOL dataset
def get_radar_from_berr(time):
    from datetime import timedelta, datetime
    import glob
    year_str = "%04d" % time.year
    month_str = "%02d" % time.month
    day_str = "%02d" % time.day
    hour_str = "%02d" % time.hour
    minute_str = "%02d" % time.minute
    second_str = "%02d" % time.second

    file_name_str = (data_path_berr_cfradial +
                     '/' + year_str + '/' + year_str + month_str +
                     day_str + '/cfrad.' + year_str + month_str +
                     day_str + '_' + hour_str + minute_str + '*.nc')
    file_name = glob.glob(file_name_str)
    print(file_name_str)
    print(file_name)
    radar = pyart.io.read(file_name[0])
    return radar


# Get a Radar object given a time period in the CPOL dataset
def get_grid_from_cpol(time):
    from datetime import timedelta, datetime
    year_str = "%04d" % time.year
    month_str = "%02d" % time.month
    day_str = "%02d" % time.day
    hour_str = "%02d" % time.hour
    minute_str = "%02d" % time.minute
    second_str = "%02d" % time.second
    file_name_str = (cpol_grid_data_path +
                     '/' + year_str + '/' + month_str + '/' + day_str +
                     '/' + 'cpol_' + year_str + month_str + day_str +
                     hour_str + minute_str + '.nc')
    print(file_name_str)
    radar = pyart.io.read_grid(file_name_str)
    return radar


def grid_radar(radar, grid_shape=(20, 301, 301), xlim=(-150000, 150000),
               ylim=(-150000, 150000), zlim=(1000, 20000), bsp=1.0,
               min_radius=750, h_factor=4.0, nb=1.5, gatefilter=False,
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
        nb=nb,
        gatefilters=[gatefilter])
    print(time.time() - bt, 'seconds to grid radar')
    return grid


def find_nearest(array, value):
    import numpy
    idx = numpy.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or
                    math.fabs(value - array[idx-1]) <
                    math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx


""" get_sounding_times
    start_year = Start year of animation
    start_month = Start month of animation
    start_day = Start day of animation
    start_hour = Start hour of animation
    end_year = End year of animation
    end_month = End month of animation
    end_day = End day of animation
    end_minute = End minute of animation
    minute_interval = Interval in minutes between scans (default is 5)
    This procedure acquires an array of sounding times between
    start_time and end_time. """


def get_sounding_times(start_year, start_month, start_day,
                       start_hour, start_minute, end_year,
                       end_month, end_day, end_hour,
                       end_minute, minute_interval=5):

    from datetime import timedelta, datetime
    start_time = datetime(start_year,
                          start_month,
                          start_day,
                          start_hour,
                          start_minute,
                          )
    end_time = datetime(end_year,
                        end_month,
                        end_day,
                        end_hour,
                        end_minute,
                        )
    deltatime = end_time - start_time

    if(deltatime.seconds > 0 or deltatime.minutes > 0):
        no_days = deltatime.days + 1
    else:
        no_days = deltatime.days

    if(start_day != end_day):
        no_days = no_days + 1

    days = np.arange(0, no_days, 1)
    print('We are about to load sounding files for ' + str(no_days) + ' days')

    # Find the list of files for each day
    cur_time = start_time

    file_list = []
    time_list = []
    for i in days:
        year_str = "%04d" % cur_time.year
        day_str = "%02d" % cur_time.day
        month_str = "%02d" % cur_time.month
        format_str = (data_path_sounding +
                      'twpsondewnpnC3.b1.' +
                      year_str +
                      month_str +
                      day_str +
                      '*custom.cdf')

        data_list = glob.glob(format_str)

        for j in range(0, len(data_list)):
            file_list.append(data_list[j])
        cur_time = cur_time + timedelta(days=1)

    # Parse all of the dates and time in the interval and add to list
    past_time = []
    for file_name in file_list:
        date_str = file_name[-26:-11]
        year_str = date_str[0:4]
        month_str = date_str[4:6]
        day_str = date_str[6:8]
        hour_str = date_str[9:11]
        minute_str = date_str[11:13]
        second_str = date_str[13:15]

        cur_time = datetime(int(year_str),
                            int(month_str),
                            int(day_str),
                            int(hour_str),
                            int(minute_str),
                            int(second_str))
        time_list.append(cur_time)

    # Sort time list and make sure time are at least xx min apart
    time_list.sort()
    time_list_sorted = deepcopy(time_list)

    time_list_final = []
    past_time = []

    for times in time_list_sorted:
        cur_time = times

        if(past_time == []):
            past_time = cur_time
        if(cur_time - past_time >= timedelta(minutes=minute_interval) and
           cur_time >= start_time and cur_time <= end_time):
            time_list_final.append(cur_time)
            past_time = cur_time

    return time_list_final


# Get a sounding object given a time period in the CPOL dataset
def get_sounding(time):
    from datetime import timedelta, datetime
    from netCDF4 import Dataset
    year_str = "%04d" % time.year
    month_str = "%02d" % time.month
    day_str = "%02d" % time.day
    hour_str = "%02d" % time.hour
    minute_str = "%02d" % time.minute
    second_str = "%02d" % time.second
    file_name_str = (data_path_sounding +
                     'twpsondewnpnC3.b1.' +
                     year_str +
                     month_str +
                     day_str +
                     '.' +
                     hour_str +
                     minute_str +
                     second_str +
                     '.custom.cdf')
    sounding = Dataset(file_name_str, mode='r')
    return sounding


# get_radar_times
#     start_year = Start year of animation
#     start_month = Start month of animation
#     start_day = Start day of animation
#     start_hour = Start hour of animation
#     end_year = End year of animation
#     end_month = End month of animation
#     end_day = End day of animation
#     end_minute = End minute of animation
#     minute_interval = Interval in minutes between scans (default is 5)
# This procedure acquires a datetime array between start_time and end_time
def get_dda_times(start_year, start_month, start_day,
                  start_hour, start_minute, end_year,
                  end_month, end_day, end_hour,
                  end_minute, minute_interval=5):
    from datetime import timedelta, datetime
    start_time = datetime(start_year, start_month, start_day,
                          start_hour, start_minute, )
    end_time = datetime(end_year, end_month, end_day, end_hour,
                        end_minute,)

    deltatime = end_time - start_time

    if(deltatime.seconds > 0 or deltatime.minute > 0):
        no_days = deltatime.days + 1
    else:
        no_days = deltatime.days

    if(start_day != end_day):
        no_days = no_days + 1

    days = np.arange(0, no_days, 1)
    print('We are about to load grid files for ' + str(no_days) + ' days')

    # Find the list of files for each day
    cur_time = start_time
    file_list = []
    time_list = []
    for i in days:
        year_str = "%04d" % cur_time.year
        day_str = "%02d" % cur_time.day
        month_str = "%02d" % cur_time.month
        format_str = (out_data_path +
                      '/ddop/cf_compliant_grid' +
                      year_str +
                      month_str +
                      day_str +
                      '*.nc')
        print('Looking for files with format ' + format_str)
        data_list = glob.glob(format_str)

        for j in range(0, len(data_list)):
            file_list.append(data_list[j])
        cur_time = cur_time + timedelta(days=1)

    # Parse all of the dates and time in the interval and add to list
    past_time = []
    for file_name in file_list:
        date_str = file_name[-15:-3]
        year_str = date_str[0:4]
        month_str = date_str[4:6]
        day_str = date_str[6:8]
        hour_str = date_str[8:10]
        minute_str = date_str[10:12]
        second_str = '00'

        cur_time = datetime(int(year_str), int(month_str), int(day_str),
                            int(hour_str), int(minute_str), 0)
        time_list.append(cur_time)

    # Sort time list and make sure time are at least xx min apart
    time_list.sort()
    time_list_sorted = deepcopy(time_list)

    time_list_final = []
    past_time = []

    for times in time_list_sorted:
        cur_time = times

        if(past_time == []):
            past_time = cur_time

        if(cur_time - past_time >= timedelta(minutes=minute_interval) and
           cur_time >= start_time and cur_time <= end_time):

            time_list_final.append(cur_time)
            past_time = cur_time

    return time_list_final


# Get a Grid object given a time period in the Multidop dataset
def get_grid_from_dda(time):
    year_str = "%04d" % time.year
    month_str = "%02d" % time.month
    day_str = "%02d" % time.day
    hour_str = "%02d" % time.hour
    minute_str = "%02d" % time.minute
    second_str = "%02d" % time.second
    file_name_str = (out_data_path +
                     '/ddop/cf_compliant_grid' +
                     year_str + month_str + day_str + hour_str +
                     minute_str + '.nc')

    radar = pyart.io.read_grid(file_name_str)
    return radar


def write_grid(time, grid):
    year_str = "%04d" % time.year
    month_str = "%02d" % time.month
    day_str = "%02d" % time.day
    hour_str = "%02d" % time.hour
    minute_str = "%02d" % time.minute
    second_str = "%02d" % time.second
    file_name_str = (out_data_path +
                     '/ddop/cf_compliant_grid' +
                     year_str + month_str + day_str + hour_str +
                     minute_str + '.nc')

    radar = grid.write(file_name_str)
    return radar
