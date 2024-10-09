'''
Python modules for making spatial plots of ELM outputs
'''
import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

__author__ = 'Eva Sinha'
__email__  = 'eva.sinha@pnnl.gov'

plt.rc('figure', titlesize=20)
plt.rc('legend', fontsize=20, title_fontsize=20)
plt.rc('axes',   labelsize=20, titlesize=20)
plt.rc('xtick',  labelsize=20)
plt.rc('ytick',  labelsize=20)
plt.rc('figure', figsize=(11, 8.5))

myDict_units = {'TOTECOSYSC':'TOTECOSYSC [PgC]',
                'TBOT':      'TBOT [K]',
                'LAND_USE_FLUX': 'LAND_USE_FLUX [PgC/year]',
                'WOOD_HARVESTC': 'WOOD_HARVESTC [PgC/year]'}

myDict_labels = {'TOTECOSYSC':'Total ecosystem carbon',
                 'LAND_USE_FLUX': 'C emissions from land use change',
                 'WOOD_HARVESTC': 'wood harvest carbon (to product pools)',
                 'TBOT':      'Atmospheric air temperature'}

# -----------------------------------------------------------
# List of variable names that we want to keep
varnames = ['TOTECOSYSC','TBOT', 'LAND_USE_FLUX', 'WOOD_HARVESTC']

mycolors=['#6A3D9A','#E31A1C','#000000','#FF7F00','#1F78B4','#000000']

for ind, var in enumerate(varnames):

  fname = 'e3sm_outputs/' + var + '_lineplot_2015-2100_global.csv'
  plot_data = pd.read_csv(fname)
  # Convert Year from float to int
  plot_data = plot_data.astype({'# Year':'int'})

  runs = list(plot_data.columns)
  runs.remove('# Year')

  for m, run in enumerate(runs):
    plt.plot(plot_data['# Year'], plot_data[run], color=mycolors[m], linewidth=0.75, marker='o')

  plt.legend(runs, bbox_to_anchor=(0.5, -0.1), ncol=4, loc='center')
  plt.ylabel(myDict_units[var])
  plt.title(myDict_labels[var])
  plt.savefig(var+'_lineplot_2015-2100_global.png', bbox_inches='tight')
  plt.clf()
