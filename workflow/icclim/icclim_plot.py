import sys
import os
import numpy as np
import xarray as xr
from scipy import stats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from icclim_levels import *

# Colormap for sptial plots
colormap = 'WhiteBlueGreenYellowRed.rgb'
rgb_arr = np.loadtxt(colormap)
rgb_arr = rgb_arr / 255.0
cmap = LinearSegmentedColormap.from_list(name=colormap, colors=rgb_arr)

# -----------------------------------------------------------
# Source: https://github.com/E3SM-Project/e3sm_diags/blob/master/e3sm_diags/driver/lat_lon_driver.py
def create_metrics(da):
    """Creates the mean, max, min, rmse, corr in a dictionary"""
    # For input None, metrics are instantiated to 999.999.
    # Apply float() to make sure the elements in metrics_dict are JSON serializable, i.e. np.float64 type is JSON serializable, but not np.float32.
    missing_value = 999.999
    metrics_dict = {}

    metrics_dict = { 
        'min': float(da.min()) if da is not None else missing_value,
        'max': float(da.max()) if da is not None else missing_value,
        #'mean': float(mean(da)) if da is not None else missing_value,# Was giving error - You have specified an invalid axis
        'mean': float(da.mean()) if da is not None else missing_value,
    }   

    return(metrics_dict)

#----------------------------------------------------------
# Calculate the T-test for the means of two independent samples of scores
def apply_ttest_ind(df, var, sel_set):

   a = df[df['Set'] == 'CONTROL'][var]
   b = df[df['Set'] == sel_set][var]

   ttest = stats.ttest_ind(a,b)

   return(ttest.pvalue)

# --------------------------------------------------------
def plot_icclim_output(out_f, index, title):
   # --- Open the saved file for plotting ---
   with xr.open_dataset('icclim_outputs/'+ out_f, decode_times=False) as ds:
      ds['time'] = xr.decode_cf(ds).time

   ds = ds[index]

   # Calculate time average
   plot_data = ds.mean(dim='time', keep_attrs=True, skipna=True)

   levels  = None
   norm    = None
   clevels = icclim_clevels[index]
   if len(clevels) > 0:
      levels = clevels
      norm = colors.BoundaryNorm(boundaries=levels, ncolors=len(levels))

   # Estimate metrics for the data
   metrics_dict = create_metrics(plot_data)
   plot_stats   = (metrics_dict['max'], metrics_dict['mean'], metrics_dict['min'])

   # Create spatial plot
   xr_plot_global(plot_data, plot_stats, \
                  cmap_col=cmap, title=title, fig_wt=8*2, fig_ht=12, \
                  levels=levels, norm=norm, fname=out_f.replace('nc', 'png'))

# --------------------------------------------------------
def plot_icclim_diff_output(out_f, out_f_r, base_case, ref_case, index, title, fname):
   # --- Open the saved file for plotting ---
   ds   = xr.open_dataset('icclim_outputs/'+ out_f, decode_times=False)
   ds_r = xr.open_dataset('icclim_outputs/'+ out_f_r, decode_times=False)

   ds['time']   = xr.decode_cf(ds).time
   ds_r['time'] = xr.decode_cf(ds_r).time

   ds   = ds[index]
   ds_r = ds_r[index]

   # Convert datset to dataarray
   ds_merge = xr.merge([ds.expand_dims(Set = [base_case]), ds_r.expand_dims(Set = [ref_case])])
   ds_merge = ds_merge.to_array()

   # Convert to pandas dataframe
   # and grouby lat and lon so that T-test can be applied on each grid
   df = ds_merge.to_dataframe(name=index).reset_index()
   df = df.dropna()
   df = df.groupby(['lat','lon'])

   # Calculate the T-test for the means of two independent samples of scores
   p_values = df.apply(apply_ttest_ind, index, sel_set=base_case)
   # Convert to xarray
   p_values = p_values.to_xarray()

   # Calculate time average
   ds = ds.mean(dim='time', keep_attrs=True, skipna=True)
   ds_r = ds_r.mean(dim='time', keep_attrs=True, skipna=True)

   # Difference
   plot_data = ds - ds_r

   # Percent difference
   plot_data_per_diff = 100*(ds - ds_r)/ds_r
   plot_data_per_diff = xr.where(np.isinf(plot_data_per_diff), 0, plot_data_per_diff)

   levels  = None
   norm    = None
   if(index == 'TXx'):
     clevels = icclim_TXx_diff_clevels
   elif(index == 'TNn'):
     clevels = icclim_TNn_diff_clevels
   elif(index == 'R20mm'):
     clevels = icclim_R20mm_diff_clevels
   else:
     clevels = icclim_diff_clevels

   if len(clevels) > 0:
      levels = clevels
      norm = colors.BoundaryNorm(boundaries=levels, ncolors=len(levels))

   # Estimate metrics for the data
   metrics_dict = create_metrics(plot_data)
   plot_stats   = (metrics_dict['max'], metrics_dict['mean'], metrics_dict['min'])

   # Create spatial plot
   xr_plot_global(plot_data, plot_stats, \
                  cmap_col='bwr', title=title, fig_wt=8*2, fig_ht=12, \
                  levels=levels, norm=norm, fname=fname,
                  stipple_data = p_values)

   levels  = None
   norm    = None
   clevels = icclim_per_diff_clevels
   if len(clevels) > 0:
      levels = clevels
      norm = colors.BoundaryNorm(boundaries=levels, ncolors=len(levels))

   # Estimate metrics for the data
   metrics_dict = create_metrics(plot_data_per_diff)
   plot_stats   = (metrics_dict['max'], metrics_dict['mean'], metrics_dict['min'])

   # Create spatial plot
   xr_plot_global(plot_data_per_diff, plot_stats, \
                  cmap_col='bwr', title=title+ ' (% diff)', fig_wt=8*2, fig_ht=12, \
                  levels=levels, norm=norm, fname=fname.replace('_diff_', '_per_diff_'),
                  stipple_data = p_values)

# -----------------------------------------------------------
def xr_plot_global(da_plot, plot_stats, cmap_col, title, fig_wt, fig_ht, levels, norm, fname, stipple_data=None):

    fig = plt.figure(figsize=(fig_wt, fig_ht))
    ax  = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection = ccrs.Robinson(central_longitude=0))

    cmap = plt.get_cmap(cmap_col)
    if(cmap_col == 'bwr'):
      cmap.set_extremes(under='darkblue', over='darkred')
    else:
      cmap.set_extremes(under='white', over='darkred')

    # spatial plot using xarray
    #fg = da_plot.plot.contourf(ax          = ax,
    fg = da_plot.plot(ax          = ax,
                               transform   = ccrs.PlateCarree(), # coordinate system of data
                               levels      = levels, 
                               #norm        = norm,
                               cmap        = cmap,
                               extend      = 'both',
                               add_colorbar= False)

    # Add stippling
    if stipple_data is not None:
      mask = stipple_data <= 0.05
      tmp = mask.values
      tmp = tmp[tmp == True]
      print('Number of grid cells with significant change')
      print(tmp.size)
      ax.contourf(stipple_data.lon, stipple_data.lat, mask,  1, hatches=['', 'xxx'], alpha=0,
                  transform=ccrs.PlateCarree())

    # Add title
    ax.set_title(title)

    # Add colorbar
    cax = fig.add_axes([ax.get_position().x1+0.01, ax.get_position().y0, 0.02, ax.get_position().height])
    cbar = plt.colorbar(fg, cax=cax)

    if levels is None:
       cbar.ax.tick_params(labelsize=15, length=0)

    else:
       maxval = np.amax(np.absolute(levels[1:-1]))
       if maxval < 0.2:
          fmt = "%5.3f"
          pad = 28
       elif maxval < 10.0:
          fmt = "%5.2f"
          pad = 40
       elif maxval < 100.0:
          fmt = "%5.0f"
          pad = 25
       elif maxval < 10000.0:
          fmt = "%5.0f"
          pad = 30
       elif maxval > 9999.0:
          fmt = "%.0f"
          pad = 50
       else:
          fmt = "%6.1f"
          pad = 40

       cbar.set_ticks(levels[1:-1])
       labels = [fmt % level for level in levels[1:-1]]
       cbar.ax.set_yticklabels(labels, ha='right')
       cbar.ax.tick_params(labelsize=15, pad=pad, length=0)

    # Min, Mean, Max
    ax.text(x=0.88, y=0.9, s='Max\nMean\nMin', ha='left', fontsize=15, transform=ax.transAxes)
    ax.text(x=1.0, y=0.9, s='%.2f\n%.2f\n%.2f'%plot_stats[0:3], ha='right', fontsize=15, transform=ax.transAxes)

    # Add additional features like coastline and oceans
    ax.coastlines(lw = 0.6, zorder=2)
    ax.add_feature(cfeature.OCEAN, facecolor='white', zorder=1)

    plt.savefig(fname, bbox_inches='tight')

    plt.close(fig=None)

    # Change directory
    os.chdir('../../workflow/icclim/')

# -----------------------------------------------------------
