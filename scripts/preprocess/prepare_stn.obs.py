#!/usr/local/anaconda/bin/python

''' This script reads in streamflow gauge data and convert to inverse routing input format '''

import datetime as dt
import pandas as pd
import argparse
import my_functions

parser = argparse.ArgumentParser()
parser.add_argument("--cfg", type=str,  help="config file for this script")
args = parser.parse_args()
cfg = my_functions.read_config(args.cfg)

start_date = dt.datetime(cfg['PARAM']['start_date'][0], \
                         cfg['PARAM']['start_date'][1], \
                         cfg['PARAM']['start_date'][2])
end_date = dt.datetime(cfg['PARAM']['end_date'][0], \
                         cfg['PARAM']['end_date'][1], \
                         cfg['PARAM']['end_date'][2])

#======================================================#
# Load station list
#======================================================#
f = open(cfg['INPUT']['list_stn_path'])
list_stn = []  # [ [station_code, lat, lon, (column)] ]
dict_stn_info = {}  # {station_code: lat, lon, (column)}
while 1:
    line = f.readline().rstrip("\n")
    if line=="":
        break
    list_stn.append(line.split()[0])
    if cfg['INPUT']['data_formst']=='USGS':
        dict_stn_info[line.split()[0]] = [float(line.split()[1]), 
                                          float(line.split()[2]), 
                                          int(line.split()[3])]
    else:
        dict_stn_info[line.split()[0]] = [float(line.split()[1]), 
                                          float(line.split()[2])]

#======================================================#
# Load data
#======================================================#
# Load data and select time range needed
dict_df_stn = {}  # a dictionary of station data
                  # {station_code: df}
for stn in list_stn:  # for each gauge station, load data
    print 'station {}'.format(stn)
    # Load data
    filename = '{}/{}'.format(cfg['INPUT']['stn_data_dir'], stn)
    if cfg['INPUT']['data_formst']=='USGS':
        column = dict_stn_info[stn][2]
        dict_df_stn[stn] = my_functions.read_USGS_data(filename, [column], ['Discharge'])
    elif cfg['INPUT']['data_formst']=='Lohmann':
        dict_df_stn[stn] = my_functions.read_Lohmann_route_daily_output(filename)

    # Select time range needed
    dict_df_stn[stn] = my_functions.select_time_range(dict_df_stn[stn], \
                                                      start_date, end_date)
    # Convert data to cfs
    if cfg['PARAM']['input_flow_unit']=='cfs':
        pass

#======================================================#
# Write basin.stn.list and basin.stn.obs
#======================================================#
# Write basin.stn.list
f = open(cfg['OUTPUT']['basin_stn_list_path'], 'w')
for i, stn in enumerate(list_stn):
    f.write('1 {} {} {} 1000 1000 1000\n'.format(i+1, dict_stn_info[stn][0], \
                                                 dict_stn_info[stn][1]+360))
f.close()

# Write basin.stn.obs
#   Put all station data into one dataframe
df_stn_all = pd.DataFrame(dict_df_stn[list_stn[0]])
df_stn_all.columns = ['stn_1']
for i in range(1,len(list_stn)):
    df_stn_all['stn_{}'.format(i+1)] = dict_df_stn[list_stn[i]]
#   Write to file
df_stn_all.to_csv(cfg['OUTPUT']['basin_stn_obs_path'], sep=' ', \
                  header=False, index=False)








