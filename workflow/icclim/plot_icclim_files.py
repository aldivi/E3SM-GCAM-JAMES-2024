import sys
import glob
import os
import datetime
import cftime
import numpy as np
import pandas as pd
import xarray as xr
import dask
import shutil

from icclim_plot import *

# --------------------------------------------------------

# --- studied period ---
base_case    = 'FULL_FDBK'
yr1       = 2071
yr2       = 2090

# --- reference period ---
ref_case  = 'CONTROL'
yr1r      = 2071
yr2r      = 2090

myDict_icclim = {
                      'TXx': {'varname'   :'TREFMXAV',
                               'temp_agg' :'year',
                               'qualifier':'',
                               'out_unit' :'degree_Celsius',
                               'title'    :'Maximum daily maximum temperature'},
                      'TNn': {'varname'   :'TREFMNAV',
                              'temp_agg' :'year',
                              'qualifier':'',
                              'out_unit' :'degree_Celsius',
                              'title'    :'Minimum daily minimum temperature'},
                      'CDD':  {'varname'  :'RAIN',
                               'temp_agg' :'year',
                               'out_unit' :'day',
                               'qualifier':'',
                               'title'    :'Maximum consecutive dry days (Precip < 1mm)'},
                      'PRCPTOT': {'varname'  :'RAIN',
                                   'temp_agg' :'year',
                                   'out_unit':'mm',
                                   'qualifier':'',
                                   'title'    :'Total precipitation during Wet Days'},
                      'RX1day':{'varname' :'RAIN',
                                'temp_agg':'year',
                                'out_unit':'mm/day',
                                'qualifier':'',
                                'title'   :'Maximum 1-day precipitation'},
                      'RX5day':{'varname' :'RAIN',
                                'temp_agg':'year',
                                'out_unit':'mm/day',
                                'qualifier':'',
                                'title'   :'Maximum 5-day precipitation'},
                      'R20mm': {'varname' :'RAIN',
                                'temp_agg':'year',
                                'out_unit':'day',
                                'qualifier':'',
                                'title'   :'Number of very heavy precipitation days (Precip >= 20mm)'},
                      'CDD':  {'varname'  :'RAIN',
                               'temp_agg' :'year',
                               'out_unit' :'day',
                               'qualifier':'',
                               'title'    :'Maximum consecutive dry days (Precip < 1mm)'}}

for i, key in enumerate(myDict_icclim):
  # ----- Estimate climate index for the reference period -----
  out_f    = key + '_' + myDict_icclim[key]['temp_agg'] + '_' + base_case + '_icclim.nc'
  out_f_r    = key + '_' + myDict_icclim[key]['temp_agg'] + '_' + ref_case + '_icclim.nc'

  fname = key + '_' + myDict_icclim[key]['temp_agg'] + '_' + base_case + '_diff_' + ref_case + '.png'
  title = myDict_icclim[key]['title'] + '\n' + base_case + ' - '+ ref_case
  plot_icclim_diff_output(out_f, out_f_r, base_case, ref_case, key, title, fname)
