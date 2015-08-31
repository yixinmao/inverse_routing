#!/usr/local/anaconda/bin/python

''' This scripts compares routed flow from inverted runoff with streamflow station observation '''

import numpy as np
import datetime as dt
import argparse
import my_functions

parser = argparse.ArgumentParser()
parser.add_argument("--cfg", type=str,  help="config file for this script")
args = parser.parse_args()
cfg = my_functions.read_config(args.cfg)

#===============================================================#
# Read data
#===============================================================#
dict_path = {}  # {station_name: [path_for_orig_obs, path_for_routed, (column in USGS data)]}
f = open(cfg['INPUT']['cmp_routed_obs_list_path'], 'r')
while 1:
    line = f.readline().rstrip("\n")
    line_split = line.split()
    if line=="":
        break
    if cfg['INPUT']['obs_format']=='USGS':
        dict_path[line_split[0]] = [line_split[1], line_split[2], int(line_split[3])]
    elif cfg['INPUT']['obs_format']=='Lohmann':
        dict_path[line_split[0]] = [line_split[1], line_split[2]]
    else:
        print 'Error: unsupported observation data format!'
        exit()
f.close()    

# Read in routed streamflow from inverted runoff
dict_Lohmann_routed = {}  # {station_name: pd.Series of daily data} [unit: cfs]
for stn in dict_path:
    # Load data
    s_Lohmann_routed = my_functions.read_Lohmann_route_daily_output(dict_path[stn][1])
    dict_Lohmann_routed[stn] = s_Lohmann_routed
    # Select full water years
    start_date_WY, end_date_WY = my_functions.find_full_water_years_within_a_range(\
                                                dict_Lohmann_routed[stn].index[0], \
                                                dict_Lohmann_routed[stn].index[-1])
    dict_Lohmann_routed[stn] = my_functions.select_time_range(dict_Lohmann_routed[stn], \
                                                              start_date_WY, \
                                                              end_date_WY)

# Read in original station obs rmat
dict_obs = {}  # {station_name: pd.Series of daily data} [unit: cfs]
for stn in dict_path:
    # Load data
    filename = dict_path[stn][0]
    if cfg['INPUT']['obs_format']=='USGS':
        column = dict_path[stn][2]
        dict_obs[stn] = my_functions.read_USGS_data(filename, [column], ['Discharge'])
    elif cfg['INPUT']['obs_format']=='Lohmann':
        dict_obs[stn] = my_functions.read_Lohmann_route_daily_output(filename)

    # Select the same range as Lohmann routed flow
    dict_obs[stn] = my_functions.select_time_range(dict_obs[stn], \
                                                   start_date_WY, \
                                                   end_date_WY)
    # Convert data to cfs
    if cfg['PARAM']['obs_flow_unit']=='cfs':
        pass

#===============================================================#
# Plot and compare
#===============================================================#
for stn in dict_path:
    # Plot daily
    fig = my_functions.plot_time_series(\
        plot_date=True, list_s_data=[dict_obs[stn], dict_Lohmann_routed[stn]], \
        list_style=['b-', 'r--'], list_label=['Obs.', 'Routed from inverse runoff'], \
        plot_start=dict_Lohmann_routed[stn].index[0], \
        plot_end=dict_Lohmann_routed[stn].index[-1], \
        xlabel=None, ylabel='Streamflow (cfs)', \
        title='Daily, {}'.format(stn), fontsize=16, legend_loc='upper right', \
        time_locator=None, time_format='%Y/%m', \
        xtick_location=None, xtick_labels=None, \
        add_info_text=False, model_info=None, stats=None, \
        show=True)
    fig.savefig('{}/cmp_obs_routed.flow_daily.{}.png'\
                        .format(cfg['OUTPUT']['output_plot_dir'], stn), \
                format='png')
    # Plot monthly
    fig = my_functions.plot_monthly_data(\
        list_s_data=[dict_obs[stn], dict_Lohmann_routed[stn]], 
        list_style=['b-', 'r--'], list_label=['Obs.', 'Routed from inverse runoff'], \
        plot_start=dict_Lohmann_routed[stn].index[0], \
        plot_end=dict_Lohmann_routed[stn].index[-1], \
        xlabel=None, ylabel='Streamflow (cfs)', \
        title='Monthly, {}'.format(stn), fontsize=16, legend_loc='upper right', \
        time_locator=None, time_format='%Y/%m', \
        add_info_text=False, model_info=None, stats=None, \
        show=False)
    fig.savefig('{}/cmp_obs_routed.flow_monthly.{}.png'\
                        .format(cfg['OUTPUT']['output_plot_dir'], stn), \
                format='png')
    # Plot seasonality
    fig = my_functions.plot_seasonality_data(\
        list_s_data=[dict_obs[stn], dict_Lohmann_routed[stn]], \
        list_style=['b-', 'r--'], list_label=['Obs.', 'Routed from inverse runoff'], \
        plot_start=1, plot_end=12, \
        xlabel=None, ylabel='Streamflow (cfs)', \
        title='Seasonality, {}, WY{}-{}'.format(stn, start_date_WY.year+1, \
                                              end_date_WY.year), \
        fontsize=16, legend_loc='upper right', \
        xtick_location=None, xtick_labels=None, \
        add_info_text=False, model_info=None, stats=None, \
        show=False)
    fig.savefig('{}/cmp_obs_routed.flow_seas.{}.png'\
                        .format(cfg['OUTPUT']['output_plot_dir'], stn), \
                format='png')






