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

def plot_time_series(plot_date, list_s_data, list_style, list_label, plot_start, plot_end, xlabel=None, ylabel=None, title=None, fontsize=16, legend_loc='lower right', time_locator=None, time_format='%Y/%m', xtick_location=None, xtick_labels=None, add_info_text=False, model_info=None, stats=None, show=False):
    ''' This function plots daily data time series

    Input:
        plot_date: True for plot_date, False for plot regular time series
        list_s_data: a list of pd.Series objects to be plotted
        list_style: a list of plotting style (e.g., ['b-', 'r--']); must be the same size as 'list_s_data'
        list_label: a list of plotting label (e.g., ['Scenario1', 'Scenario2']); must be the same size as 'list_s_data'
        xlabel: [str]
        ylabel: [str]
        title: [str]
        fontsize: for xlabe, ylabel and title [int]
        legend_loc: [str]
        plot_start, plot_end: if plot_date=True, [dt.datetime]; if plot_date=False, [float/int]
        time_locator: time locator on the plot; 'year' for year; 'month' for month. e.g., ('month', 3) for plot one tick every 3 months [tuple]
        time_format: [str]
        xtick_location: a list of xtick locations [list of float/int]
        xtick_labels: a list of xtick labels [list of str]; must be the same length as 'xtick_locations'
        add_info_text: True for adding info text at the bottom of the plot
        model_info, stats: descriptions added in the info text [str]
        show: True for showing the plot

    Require:
        plot_date_format
        add_info_text_to_plot(fig, ax, model_info, stats)
        plot_format
    '''

    import matplotlib.pyplot as plt

    # Check if list_s_data, list_style and list_label have the same length
    if len(list_s_data) !=len(list_style) or len(list_s_data)!=len(list_label):
        print 'Input list lengths are not the same!'
        exit()

    fig = plt.figure(figsize=(12,8))
    ax = plt.axes()
    for i in range(len(list_s_data)):
        if plot_date==True:  # if plot date
            plt.plot_date(list_s_data[i].index, list_s_data[i], list_style[i], label=list_label[i])
        else:  # if plot regular time series
            plt.plot(list_s_data[i].index, list_s_data[i], list_style[i], label=list_label[i])
    if xlabel:
        plt.xlabel(xlabel, fontsize=fontsize)
    if ylabel:
        plt.ylabel(ylabel, fontsize=fontsize)
    if title:
        plt.title(title, fontsize=fontsize)
    # format plot
    leg = plt.legend(loc=legend_loc, frameon=True)
    leg.get_frame().set_alpha(0)
    if plot_date==True:  # if plot date
        plot_date_format(ax, time_range=(plot_start, plot_end), locator=time_locator, time_format='%Y/%m')
    else:  # if plot regular time series
        plt.xlim([plot_start, plot_end])
        if xtick_location:
            plot_format(ax, xtick_location=xtick_location, xtick_labels=xtick_labels)
    # add info text
    if add_info_text==True:
        add_info_text_to_plot(fig, ax, model_info, stats)

    if show==True:
        plt.show()

    return fig

#==============================================================
#==============================================================

def plot_date_format(ax, time_range=None, locator=None, time_format=None):
    ''' This function formats plots by plt.plot_date

    Input:
        ax: plotting axis
        time_range: a tuple of two datetime objects indicating xlim. e.g., (dt.date(1991,1,1), dt.date(1992,12,31))
        locator: time locator on the plot; 'year' for year; 'month' for month. e.g., ('month', 3) for plot one tick every 3 months
        time_format: a string of time format, e.g. '%Y/%m'
    '''

    import matplotlib.pyplot as plt
    import datetime as dt
    from matplotlib.dates import YearLocator, MonthLocator, DateFormatter

    # Plot time range
    if time_range!=None:
        plt.xlim(time_range[0], time_range[1])

    # Set time locator (interval)
    if locator!=None:
        if locator[0]=='year':
            ax.xaxis.set_major_locator(YearLocator(locator[1]))
        elif locator[0]=='month':
            ax.xaxis.set_major_locator(MonthLocator(interval=locator[1]))

    # Set time ticks format
    if time_format!=None:
        ax.xaxis.set_major_formatter(DateFormatter(time_format))

    return ax

#==============================================================
#==============================================================

def add_info_text_to_plot(fig, ax, model_info, stats, fontsize=14, bottom=0.3, text_location=-0.1):
    ''' This function adds info text to the bottom of a plot
        The text will include:
            Author;
            Plotting date;
            Model info (taking from input [str]);
            Stats (taking from input [str])
            bottom: the extent of adjusting the original plot; the higher the 'bottom', the more space it would be left for the texts
            text_location: the location of the text; the more negative, the lower the textG
'''

    import matplotlib.pyplot as plt
    import datetime as dt

    # adjust figure to leave some space at the bottom
    fig.subplots_adjust(bottom=bottom)

    # determine text content
    author = 'Yixin'
    today = dt.date.today()
    plot_date = today.strftime("%Y-%m-%d")
    text_to_add = 'Author: %s\nDate plotted: %s\nModel info: %s\nStats: %s\n' %(author, plot_date, model_info, stats)

    # add text
    plt.text(0, text_location, text_to_add, horizontalalignment='left',\
            verticalalignment='top', transform=ax.transAxes, fontsize=fontsize)

    return fig, ax

#==============================================================
#==============================================================

def plot_format(ax, xtick_location=None, xtick_labels=None):
    '''This function formats plots by plt.plot

    Input:
        xtick_location: e.g. [1, 2, 3]
        xtick_labels: e.g. ['one', 'two', 'three']
    '''

    import matplotlib.pyplot as plt

    ax.set_xticks(xtick_location)
    ax.set_xticklabels(xtick_labels)

    return ax

#==============================================================
#==============================================================

def plot_monthly_data(list_s_data, list_style, list_label, plot_start, plot_end, xlabel=None, ylabel=None, title=None, fontsize=16, legend_loc='lower right', time_locator=None, time_format='%Y/%m', add_info_text=False, model_info=None, stats=None, show=False):
    ''' This function plots monthly mean data time series

    Require:
        plot_date_format
        add_info_text_to_plot(fig, ax, model_info, stats)
        plot_time_series
        calc_monthly_data
    '''

    # Check if list_s_data, list_style and list_label have the same length
    if len(list_s_data) !=len(list_style) or len(list_s_data)!=len(list_label):
        print 'Input list lengths are not the same!'
        exit()

    # Calculate monthly mean data
    list_s_month = []   # list of new monthly mean data in pd.Series type
    for i in range(len(list_s_data)):
        s_month = calc_monthly_data(list_s_data[i])
        list_s_month.append(s_month)

    # plot
    fig = plot_time_series(True, list_s_month, list_style, list_label, plot_start, plot_end, xlabel, ylabel, title, fontsize, legend_loc, time_locator, time_format, add_info_text=add_info_text, model_info=model_info, stats=stats, show=show)

    return fig

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

def generate_xmask_for_route(flowdir_file):
    ''' This function generates xmask (i.e., flow distance) data using haversine formula
    Input:
        Flow direction file path, in the format of 1-8 and 9 for outlet
    Return:
        A np.array matrix of flow distance [unit: m], -1 for inactive grid cells

    Require:
        read_GIS_ascii_header        
    '''

    import numpy as np

    r_earth = 6371.0072 * 1000  # earth radius [unit: m]

    #=== Read in flow direction file ===#
    # Read header
    ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value = \
            read_GIS_ascii_header(flowdir_file)
    # Read flow direction
    fdir_all = np.loadtxt(flowdir_file, dtype=int, skiprows=6)

    #=== Loop over each grid cell in column order ===#
    flow_dist_grid = np.ones([nrows, ncols]) * -1
    lat_max = yllcorner + nrows*cellsize - cellsize/2.0  # Northmost grid cell lat
    for j in range(ncols):
    # Grid cell lon
        grid_lon = xllcorner + cellsize/2.0 + j*cellsize
        for i in range(nrows):
            # Grid cell lat
            grid_lat = lat_max - i*cellsize
            # Calculate flow distance, if active cell
            if fdir_all[i][j]!=int(NODATA_value):  # if active cell
                # Get flow direction
                fdir = fdir_all[i][j]
                # Determine lat and lon of 1st order downstream grid cell
                if fdir==1 or fdir==2 or fdir==8:
                    ds1_lat = grid_lat + cellsize
                elif fdir==4 or fdir==5 or fdir==6:
                    ds1_lat = grid_lat - cellsize
                else:
                    ds1_lat = grid_lat
                if fdir==2 or fdir==3 or fdir==4:
                    ds1_lon = grid_lon + cellsize
                elif fdir==6 or fdir==7 or fdir==8:
                    ds1_lon = grid_lon - cellsize
                else:
                    ds1_lon = grid_lon
                # Calculate flow distance to the downstream grid cell
                hslat = (1 - np.cos((grid_lat-ds1_lat)/180.0*np.pi) ) / 2.0
                hslon = (1 - np.cos((grid_lon-ds1_lon)/180.0*np.pi) ) / 2.0
                flow_dist = 2 * r_earth * np.arcsin(np.sqrt(hslat + \
                                                    np.cos(grid_lat/180.0*np.pi) \
                                                    * np.cos(ds1_lat/180.0*np.pi) * hslon))
                # Save flow distance to array
                flow_dist_grid[i,j] = flow_dist

    return flow_dist_grid

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

def plot_seasonality_data(list_s_data, list_style, list_label, plot_start, plot_end, xlabel=None, ylabel=None, title=None, fontsize=16, legend_loc='lower right', xtick_location=None, xtick_labels=None, add_info_text=False, model_info=None, stats=None, show=False):
    ''' This function plots seasonality data time series (12 month's mean)

    Require:
        plot_date_format
        add_info_text_to_plot(fig, ax, model_info, stats)
        plot_time_series
        calc_ts_stats_by_group
        plot_format
    '''

    # Check if list_s_data, list_style and list_label have the same length
    if len(list_s_data) !=len(list_style) or len(list_s_data)!=len(list_label):
        print 'Input list lengths are not the same!'
        exit()

    # Calculate monthly mean data
    list_s_seas = []   # list of new monthly mean data in pd.Series type
    for i in range(len(list_s_data)):
        s_seas = calc_ts_stats_by_group(list_s_data[i], 'month', 'mean') # index is 1-12 (month)
        list_s_seas.append(s_seas)

    # plot
    fig = plot_time_series(False, list_s_seas, list_style, list_label, plot_start, plot_end, xlabel, ylabel, title, fontsize, legend_loc, xtick_location=xtick_location, xtick_labels=xtick_labels, add_info_text=add_info_text, model_info=model_info, stats=stats, show=show)

    return fig

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




