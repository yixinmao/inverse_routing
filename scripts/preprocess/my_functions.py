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

#==============================================================
#==============================================================

def read_USGS_data(file, columns, names):
    '''This function reads USGS streamflow from the directly downloaded format (date are in the 3rd column)

    Input:
        file: directly downloaded streamflow file path [str]
        columns: a list of data colomn numbers, starting from 1. E.g., if the USGS original data has three variables: max_flow, min_flow, mean_flow, and the desired variable is mean_flow, then columns = [3]
        names: a list of data column names. E.g., ['mean_flow']; must the same length as columns

    Return:
        a pd.DataFrame object with time as index and data columns (NaN for missing data points)

    Note: returned data and flow might not be continuous if there is missing data!!!

    '''

    import numpy as np
    import datetime as dt
    import pandas as pd

    ndata = len(columns)
    if ndata != len(names):  # check input validity
        print "Error: input arguments 'columns' and 'names' must have same length!"
        exit()

    f = open(file, 'r')
    date_array = []
    data = []
    for i in range(ndata):
        data.append([])
    while 1:
        line = f.readline().rstrip("\n")  # read in one line
        if line=="":
            break
        line_split = line.split('\t')
        if line_split[0]=='USGS':  # if data line
            date_string = line_split[2]  # read in date string
            date = dt.datetime.strptime(date_string, "%Y-%m-%d")  # convert date to dt object
            date_array.append(date)

            for i in range(ndata):  # for each desired data variable
                col = columns[i]
                if line_split[3+(col-1)*2] == '':  # if data is missing
                    value = np.nan
                elif line_split[3+(col-1)*2] == 'Ice':  # if data is 'Ice'
                    value = np.nan
                else:  # if data is not missing
                    value = float(line_split[3+(col-1)*2])
                data[i].append(value)

    data = np.asarray(data).transpose()
    df = pd.DataFrame(data, index=date_array, columns=names)
    return df

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

def read_GIS_ascii_header(file):
    ''' This function reads GIS ascii file header 
    Input:
        file: ascii file path; the first 6 lines are the header
    '''

    f = open(file, 'r')
    # ncols
    line = f.readline().rstrip("\n")
    if line.split()[0]!='ncols':
        print 'Error: {} - flow direction file header variable name \
               unsupported!'.format(line.split()[0])
        exit()
    ncols = int(line.split()[1])
    # nrows
    line = f.readline().rstrip("\n")
    if line.split()[0]!='nrows':
        print 'Error: {} - flow direction file header variable name \
               unsupported!'.format(line.split()[0])
        exit()
    nrows = int(line.split()[1])
    # xllcorner
    line = f.readline().rstrip("\n")
    if line.split()[0]!='xllcorner':
        print 'Error: {} - flow direction file header variable name \
               unsupported!'.format(line.split()[0])
        exit()
    xllcorner = float(line.split()[1])
    # yllcorner
    line = f.readline().rstrip("\n")
    if line.split()[0]!='yllcorner':
        print 'Error: {} - flow direction file header variable name \
               unsupported!'.format(line.split()[0])
        exit()
    yllcorner = float(line.split()[1])
    # cellsize
    line = f.readline().rstrip("\n")
    if line.split()[0]!='cellsize':
        print 'Error: {} - flow direction file header variable name \
               unsupported!'.format(line.split()[0])
        exit()
    cellsize = float(line.split()[1])
    # NODATA_value
    line = f.readline().rstrip("\n")
    if line.split()[0]!='NODATA_value':
        print 'Error: {} - flow direction file header variable name \
               unsupported!'.format(line.split()[0])
        exit()
    NODATA_value = float(line.split()[1])
    
    return ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value

#==============================================================
#==============================================================

def read_VIC_output(file, data_columns, data_names, header=False, date_col=3):
    ''' This function reads VIC-output-formatted data

    Input:
        file: VIC output file path
        header: True for with header lines (6 lines); False for without header lines
        date_col: 3 if the first three columns are datetime (year, month, day);
                  4 if the first four columns are datetime (year, month, day, hour)
        data_columns: list of columns of data wanted, index starting from 1
                      e.g., [4, 5] or [4]
        data_names: list of names of data, same length as data_columns
                    e.g., ['runoff', 'baseflow']

    Return:
        a dataframe
    '''

    import pandas as pd
    import datetime as dt

    if header==True:
        skiprows = 6
    else:
        skiprows = 0


    # load data
    if date_col==3:
        parse = lambda x: dt.datetime.strptime(x, '%Y\t%m\t%d')
        df = pd.read_csv(file, delim_whitespace=True, parse_dates=[[0,1,2]], \
                         index_col=0, date_parser=parse, header=None, \
                         skiprows=skiprows, usecols=[0,1,2]+[i-1 for i in data_columns])
    elif date_col==4:
        parse = lambda x: dt.datetime.strptime(x, '%Y\t%m\t%d\t%h')
        df = pd.read_csv(file, delim_whitespace=True, parse_dates=[[0,1,2,3]], \
                         index_col=0, date_parser=parse, header=None, \
                         skiprows=skiprows, usecols=[0,1,2]+[i-1 for i in data_columns])
    else:
        print 'Error: unsupported date column number: {}!'.format(date_col)
        exit()

    # Rename columns
    if len(data_columns)!=len(data_names):
        print 'Error: length of data columns and names not equal!'
        exit()
    for i in range(len(data_columns)):
        df = df.rename(columns={data_columns[i]-1:data_names[i]})
    
    return df

#==============================================================
#==============================================================

def read_Lohmann_route_daily_output(path):
    ''' This function reads Lohmann routing model daily output

    Input: daily output file path
    Return: a pd.Series object with datetime as index and flow[cfs] as data

    '''

    import pandas as pd
    import datetime as dt

    parse = lambda x: dt.datetime.strptime(x, '%Y %m %d')

    # load data
    df = pd.read_csv(path, delim_whitespace=True, parse_dates=[[0,1,2]], index_col=0, date_parser=parse, header=None)
    df = df.rename(columns={3:'flow'})
    # convert data to pd.Series
    s = df.ix[:,0]

    return s



