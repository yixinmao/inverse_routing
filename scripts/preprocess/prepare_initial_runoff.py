#!/usr/local/anaconda/bin/python

''' This script reads in VIC output and convert to initial runoff field guess for inverse routing input '''

import numpy as np
import pandas as pd
import datetime as dt
import argparse
import xray
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
# Load and process VIC output netCDF file
#=====================================================#
# Load data
print 'Loading data...'
ds = xray.open_dataset(cfg['INPUT']['vic_output_nc'])
# Select time period needed
print 'Selecting time range needed...'
ds = ds.sel(time=slice(start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d')))
# Sum all runoff variables
print 'Summing up all runoff variables needed...'
if type(cfg['INPUT']['vic_runoff_variables']) is str:   # if only one variable
    da = ds[cfg['INPUT']['vic_runoff_variables']]
elif type(cfg['INPUT']['vic_runoff_variables']) is list:   # if multiple variable, sum up
    for i, var in enumerate(cfg['INPUT']['vic_runoff_variables']):
        if i==0:  # for the first variable, create da
            da = ds[var]
        else:  # for the following variables
            da = da + ds[var]

#=====================================================#
# Select grid cells in the flow direction area
#=====================================================#
# Prepare a pd.Series with all -1 (for inactive grid cells), and count # days
s_inactive = pd.Series(-1, index=pd.date_range(start_date, end_date, freq='24H'))
array_inactive = s_inactive.values
nday = len(s_inactive)
ncell = ncols * nrows
# initialize nd.array
array_vic_output = np.empty([nday, ncell])
# Loop over each grid cell in column order
lat_max = yllcorner + nrows*cellsize - cellsize/2.0  # Northmost grid cell lat
count = 0
for j in range(ncols):
    # Grid cell lon
    lon = xllcorner + cellsize/2.0 + j*cellsize
    for i in range(nrows):
        print 'Loading row {}, col {}...'.format(i+1, j+1)
        # Grid cell lat
        lat = lat_max - i*cellsize
        # Extract data at this grid cell, if this is an active cell
        if fdir[i][j]!=int(NODATA_value):  # if active cell
            array_vic_output[:, count] = da.sel(lat=lat, lon=lon).values
        else:  # if inactive cell
            array_vic_output[:, count] = array_inactive

        count = count + 1

#=====================================================#
# Write runoff field to file
#=====================================================#
np.savetxt(cfg['OUTPUT']['basin_runoff_path'], array_vic_output, fmt='%.4f')




