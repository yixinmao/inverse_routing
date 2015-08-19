#!/usr/local/anaconda/bin/python

# -------------------------------------------------------------------- #
def read_config(config_file, default_config=None):
    """
    This function is from tonic (author: Joe Hamman)
    Return a dictionary with subdictionaries of all configFile options/values
    """

    from netCDF4 import Dataset
    try:
        from cyordereddict import OrderedDict
    except:
        from collections import OrderedDict
    try:
        from configparser import SafeConfigParser
    except:
        from ConfigParser import SafeConfigParser
    import configobj

    config = SafeConfigParser()
    config.optionxform = str
    config.read(config_file)
    sections = config.sections()
    dict1 = OrderedDict()
    for section in sections:
        options = config.options(section)
        dict2 = OrderedDict()
        for option in options:
            dict2[option] = config_type(config.get(section, option))
        dict1[section] = dict2

    if default_config is not None:
        for name, section in dict1.items():
            if name in default_config.keys():
                for option, key in default_config[name].items():
                    if option not in section.keys():
                        dict1[name][option] = key

    return dict1
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def config_type(value):
    """
    This function is originally from tonic (author: Joe Hamman); modified so that '\' is considered as an escapor. For example, '\,' can be used for strings with ','. e.g., Historical\, 1980s  will be recognized as one complete string
    Parse the type of the configuration file option.
    First see the value is a bool, then try float, finally return a string.
    """

    import cStringIO
    import csv

    val_list = [x.strip() for x in csv.reader(cStringIO.StringIO(value), delimiter=',', escapechar='\\').next()]
    if len(val_list) == 1:
        value = val_list[0]
        if value in ['true', 'True', 'TRUE', 'T']:
            return True
        elif value in ['false', 'False', 'FALSE', 'F']:
            return False
        elif value in ['none', 'None', 'NONE', '']:
            return None
        else:
            try:
                return int(value)
            except:
                pass
            try:
                return float(value)
            except:
                return value
    else:
        try:
            return list(map(int, val_list))
        except:
            pass
        try:
            return list(map(float, val_list))
        except:
            return val_list
# -------------------------------------------------------------------- #

def read_inverse_route_output(output_dir, smooth_window, skip_steps, start_date, time_step, latlon_precision):
    ''' This function reads in inverse routing output 
    Input:
        output_dir: inverse routing output directory (e.g., './output/basin')
        smooth_window: length of smooth window
        skip_steps: number of time steps to skip
        start_date: start date of inverse routing input (before skipping days) [dt.datetime]
        time_step: inverse routing time step [unit: hour]
        latlon_precision: number of figures after decimal point for lat and lon

    Return:
        a dictionary; keys: 'lat_lon'; element: Series of total runoff
    '''

    import numpy as np
    import pandas as pd
    import datetime as dt
    
    filename = '{}/data_all_day_{:d}_skip{:d}'.format(output_dir, smooth_window, skip_steps)
    dict_s = {} # a dictionary; keys: 'lat_lon'; element: Series of total runoff
    f = open(filename, 'r')
    while 1:
        line = f.readline().rstrip("\n")
        line_split = line.split('\t')
        if line=="":
            break
        
        # Get lat lon info
        lat_lon = '{:.{}f}_{:.{}f}'.format(float(line_split[1]), latlon_precision, \
                                           float(line_split[0]), latlon_precision)
        # Convert data to pd.Series
        total_runoff = np.asarray([float(i) for i in line_split[2:]]) # data
        end_date = start_date + dt.timedelta(hours=time_step)*(len(total_runoff)-1)
        index = pd.date_range(start_date, end_date, freq='{:d}H'.format(time_step))  # index
        dict_s[lat_lon] = pd.Series(total_runoff, index=index)
    f.close()

    return dict_s

#==============================================================
#==============================================================

def find_full_water_years_within_a_range(dt1, dt2):
    ''' This function determines the start and end date of full water years within a time range

    Input:
        dt1: time range starting time [dt.datetime]
        dt2: time range ending time [dt.datetime]

    Return:
        start and end date of full water years
    '''

    import datetime as dt

    if dt1.month <= 9:  # if dt1 is before Oct, start from Oct 1 this year
        start_date_WY = dt.datetime(dt1.year, 10, 1)
    elif dt1.month==10 and dt1.day==1:  # if dt1 is on Oct 1, start from this date
        start_date_WY = dt.datetime(dt1.year, 10, 1)
    else:  # if dt1 is after Oct 1, start from Oct 1 next year
        start_date_WY = dt.datetime(dt1.year+1, 10, 1)

    if dt2.month >=10:  # if dt2 is Oct or after, end at Sep 30 this year
        end_date_WY = dt.datetime(dt2.year, 9, 30)
    elif dt2.month==9 and dt2.day==30:  # if dt2 is on Sep 30, end at this date
        end_date_WY = dt.datetime(dt2.year, 9, 30)
    else:  # if dt2 is before Sep 30, end at Sep 30 last year
        end_date_WY = dt.datetime(dt2.year-1, 9, 30)

    if (end_date_WY-start_date_WY).days > 0:  # if at least one full water year
        return start_date_WY, end_date_WY
    else: # else, return -1
        return -1

#==============================================================
#==============================================================

def select_time_range(data, start_datetime, end_datetime):
    ''' This function selects out the part of data within a time range

    Input:
        data: [dataframe/Series] data with index of datetime
        start_datetime: [dt.datetime] start time
        end_datetime: [dt.datetime] end time

    Return:
        Selected data (same object type as input)
    '''

    import datetime as dt

    start = data.index.searchsorted(start_datetime)
    end = data.index.searchsorted(end_datetime)

    data_selected = data.ix[start:end+1]

    return data_selected

#==============================================================
#==============================================================

def calc_monthly_data(data):
    '''This function calculates monthly mean values

    Input: [DataFrame/Series] with index of time
    Return: a [DataFrame/Series] object, with monthly mean values (the same units as input data)
    '''

    import pandas as pd
    data_mon = data.resample("M", how='mean')
    return data_mon

#==============================================================
#==============================================================

def wateryear(calendar_date):
    if calendar_date.month >= 10:
        return calendar_date.year+1
    return calendar_date.year

def calc_ts_stats_by_group(data, by, stat):
    '''This function calculates statistics of time series data grouped by year, month, etc

    Input:
        df: a [pd.DataFrame/Series] object, with index of time
        by: string of group by, (select from 'year' or 'month' or 'WY')
        stat: statistics to be calculated, (select from 'mean')
        (e.g., if want to calculate monthly mean seasonality (12 values), by='month' and stat='mean')

    Return:
        A [dateframe/Series] object, with group as index (e.g. 1-12 for 'month')

    Require:
        wateryear
        find_full_water_years_within_a_range(dt1, dt2)
        select_time_range(data, start_datetime, end_datetime)
    '''

    import pandas as pd

    if by=='year':
        if stat=='mean':
            data_result = data.groupby(lambda x:x.year).mean()
    elif by=='month':
        if stat=='mean':
            data_result = data.groupby(lambda x:x.month).mean()
    elif by=='WY':  # group by water year
        # first, secelect out full water years
        start_date, end_date = find_full_water_years_within_a_range(data.index[0], data.index[-1])
        data_WY = select_time_range(data, start_date, end_date)
        # then, group by water year
        if stat=='mean':
            data_result = data_WY.groupby(lambda x:wateryear(x)).mean()

    return data_result





