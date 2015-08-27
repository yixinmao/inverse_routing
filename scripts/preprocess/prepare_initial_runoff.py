#!/usr/local/anaconda/bin/python

''' This script reads in VIC output and convert to initial runoff field guess for inverse routing input '''

import numpy as np
import pandas as pd
import datetime as dt
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

#=====================================================#
# Read in flow direction file
#=====================================================#
#----- Read header file -----#
ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value = \
    my_functions.read_GIS_ascii_header(cfg['INPUT']['fdir_header_path'])
#----- Read flow direction file -----#
fdir = np.loadtxt(cfg['INPUT']['fdir_path'], dtype=int)

#=====================================================#
# Read in VIC output files
#=====================================================#
# initialize dataframe
df_vic_output = pd.DataFrame()
# Prepare a pd.Series with all -1 (for inactive grid cells)
s_inactive = pd.Series(-1, index=pd.date_range(start_date, end_date, freq='24H'))
# Loop over each grid cell in column order
lat_max = yllcorner + nrows*cellsize - cellsize/2.0  # Northmost grid cell lat
for j in range(ncols):
    # Grid cell lon
    lon = xllcorner + cellsize/2.0 + j*cellsize
    for i in range(nrows):
        print 'Loading row {}, col {}...'.format(i+1, j+1)
        # Grid cell lat
        lat = lat_max - i*cellsize
        # Load VIC output file, if this is an active cell
        if fdir[i][j]!=int(NODATA_value):  # if active cell
            vic_output_file = '{}/{}_{:.{}f}_{:.{}f}'\
                                .format(cfg['INPUT']['vic_output_dir'], \
                                        cfg['INPUT']['vic_output_file_prefix'], \
                                        lat, cfg['INPUT']['vic_output_precise'], \
                                        lon, cfg['INPUT']['vic_output_precise'])
            # Load file
            df = my_functions.read_VIC_output(vic_output_file, data_columns=[6,7], \
                                         data_names=['runoff','baseflow'], \
                                         header=True, date_col=3)
            # Calculate total runoff
            s_total_runoff = df['runoff'] + df['baseflow']
            # Select time range interested
            s_total_runoff = my_functions.select_time_range(s_total_runoff, \
                                                            start_date, \
                                                            end_date)
            # Add to final dataframe
            df_vic_output['row{}_col{}'.format(i+1,j+1)] = s_total_runoff
                                        
        else:  # if inactive cell
            df_vic_output['row{}_col{}'.format(i+1,j+1)] = s_inactive

#=====================================================#
# Write runoff field to file
#=====================================================#
df_vic_output.to_csv(cfg['OUTPUT']['basin_runoff_path'], sep=' ', \
                  header=False, index=False)




