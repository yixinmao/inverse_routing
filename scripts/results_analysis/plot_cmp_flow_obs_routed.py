#!/usr/local/anaconda/bin/python

''' This scripts compares routed flow from inverted runoff with streamflow station observation '''

import numpy as np
import datetime as dt
import pandas as pd
import argparse
import my_functions

parser = argparse.ArgumentParser()
parser.add_argument("--cfg", type=str,  help="config file for this script")
args = parser.parse_args()
cfg = my_functions.read_config(args.cfg)

#===============================================================#
# Read station list
#===============================================================#
obs_filename_list = []
lohm_filename_list = []

f = open(cfg['INPUT']['site_list_path'])
while 1:
    line = f.readline().rstrip("\n")
    if line=="":
        break

    obs_filename_list.append(line.split()[0])
    lohm_filename_list.append('{}.day'.format(line.split()[1]))

f.close()

#===============================================================#
# Plot and compare
#===============================================================#
for i in range(len(obs_filename_list)):
    obs_filename = obs_filename_list[i]
    lohm_filename = lohm_filename_list[i]
    stn = obs_filename
    print 'Plotting {}...'.format(stn)

    # Load Lohmann routed flow
    s_Lohmann_routed = my_functions.read_Lohmann_route_daily_output('{}/{}'.format(cfg['INPUT']['Lohmann_output_dir'], lohm_filename))
    # Load obs. flow
    if cfg['INPUT']['obs_format']=='Lohmann':
        s_obs = my_functions.read_Lohmann_route_daily_output('{}/{}'.format(cfg['INPUT']['obs_dir'], obs_filename))

    # Select full water years from Lohmann routed flow
    start_date_WY, end_date_WY = my_functions.find_full_water_years_within_a_range(\
                                                s_Lohmann_routed.index[0], \
                                                s_Lohmann_routed.index[-1])
    s_Lohmann_routed = s_Lohmann_routed.truncate(before=start_date_WY, after=end_date_WY)
    s_obs = s_obs.truncate(before=start_date_WY, after=end_date_WY)

    # Plot daily
    fig = my_functions.plot_time_series(\
        plot_date=True, list_s_data=[s_obs, s_Lohmann_routed], \
        list_style=['b-', 'r--'], list_label=['Obs.', 'Routed from inverse runoff'], \
        plot_start=dt.datetime(1992,10,1),  # dict_Lohmann_routed[stn].index[0], \
        plot_end=dt.datetime(1993,9,30),   # dict_Lohmann_routed[stn].index[-1], \
        xlabel=None, ylabel='Streamflow (cfs)', \
        title='Daily, {}'.format(stn), fontsize=16, legend_loc='upper right', \
        time_locator=None, time_format='%Y/%m', \
        xtick_location=None, xtick_labels=None, \
        add_info_text=False, model_info=None, stats=None, \
        show=False)
    fig.savefig('{}/cmp_obs_routed.flow_daily.{}.png'\
                        .format(cfg['OUTPUT']['output_plot_dir'], stn), \
                format='png')
    # Plot monthly
    fig = my_functions.plot_monthly_data(\
        list_s_data=[s_obs, s_Lohmann_routed], 
        list_style=['b-', 'r--'], list_label=['Obs.', 'Routed from inverse runoff'], \
        plot_start=dt.datetime(1993,10,1),   # dict_Lohmann_routed[stn].index[0], \
        plot_end=dt.datetime(1994,9,30),   # dict_Lohmann_routed[stn].index[-1], \
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
        list_s_data=[s_obs, s_Lohmann_routed], \
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






