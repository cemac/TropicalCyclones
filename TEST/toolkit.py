# -*- coding: utf-8 -*-
"""toolkit

.. module:: toolkit
    :platform: Unix
    :synopis: Standalone tools

.. moduleauthor: John Ashcroft, CEMAC (UoL) February 2019.

.. description: This module was developed by CEMAC as part of the WCSSP
                Project. Intial script improvements

   :copyright: © 2019 University of Leeds.
   :license: BSD-2 Clause.

Example:
    To use::

Memebers:

.. CEMAC_TropicalCyclones:
   https://github.com/cemac/TropicalCyclones
"""
import numpy as np
import iris
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib
from matplotlib.colors import from_levels_and_colors
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.cm as mpl_cm


def find_subplot_dims(n):
    n_cols = np.ceil(n**0.5).astype(int)
    n_rows = np.ceil(1. * n / n_cols).astype(int)
    return n_rows, n_cols


def plot_wspeed(ax, ws):
    levels = np.arange(11) * 5
    cmap = mpl_cm.get_cmap('plasma')
    x = ws.coord('longitude').points
    y = ws.coord('latitude').points
    X, Y = np.meshgrid(x, y)

    wspeed = ax.contourf(X, Y, ws.data, cmap=cmap,
                         levels=levels, extend='both')
    return wspeed


def plot_winds(ax, u, v, mod):
    if mod == '4p4':
        scale = 300.
        step = 12
    else:
        scale = 400.
        step = 8  # Need to trial and error these.

    x = u.coord('longitude').points
    y = u.coord('latitude').points
    X, Y = np.meshgrid(x, y)

    arrows = ax.quiver(X[::step, ::step], Y[::step, ::step],
                       u.data[::step, ::step],
                       v.data[::step, ::step], angles='xy', scale=scale,
                       transform=ccrs.PlateCarree())
    if u.coord('ensemble_member').points[0] == 0:
        ax.quiverkey(arrows, 0.85, 0.1, 20, '20ms$^{-1}$',
                     coordinates='figure', fontproperties={'size': 15})


def map_formatter(ax, tick_base_x=15.0, tick_base_y=15.0, labelsize=20,
                  top_label=False, bottom_label=True, right_label=False,
                  left_label=True, central_longitude=0.0, res='10m'):
    '''
    Adds gridlines, countries and lat/lon labels to the plot.
    '''
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.75, color='k', linestyle=':')
    gl.xlabels_top = top_label
    gl.ylabels_right = right_label
    gl.ylabels_left = left_label
    gl.xlabels_bottom = bottom_label
    gl.xlocator = mticker.MultipleLocator(base=tick_base_x)
    gl.ylocator = mticker.MultipleLocator(base=tick_base_y)
    ax.coastlines(resolution=res, color='k', linewidth=1)
    gl.xlabel_style = {'size': labelsize}
    gl.ylabel_style = {'size': labelsize}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER


def find_centre_3h(md, TT, em, v_day, v_time, mod):
    '''
    Loads the track data for a specific time. This function is tailored to my
    own data, will have to be rewritten depending on how the users track data
    is stored. If you want just a single domain simply create a function that
    returns the centre of this domain.
    '''
    if mod == '4p4':
        filepath = "/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver3hr_{3}_{0:04d}_{1:02d}Z_en{2:02d}.npz".format(
            md, TT, em, mod)
    elif mod == 'glm':
        filepath = "/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver3hr_newinterp_{3}_{0:04d}_{1:02d}Z_en{2:02d}.npz".format(
            md, TT, em, mod)
    track_data = np.load(filepath)
    lats = track_data['lats']
    lons = track_data['lons']
    days = track_data['vt_days']
    times = track_data['vt_times']
    index1 = np.where(days == v_day)
    lats = lats[index1]
    lons = lons[index1]
    times = times[index1]
    index2 = np.where(times == v_time)
    cenlat = lats[index2]
    cenlon = lons[index2]
    return cenlat[0], cenlon[0]


def box_constraint(minlat, maxlat, minlon, maxlon):
    # Create constraint to extract data from cube over a certain region
    longitude_constraint1 = iris.Constraint(
        longitude=lambda cell: cell > minlon)
    longitude_constraint2 = iris.Constraint(
        longitude=lambda cell: cell < maxlon)
    latitude_constraint1 = iris.Constraint(latitude=lambda cell: cell > minlat)
    latitude_constraint2 = iris.Constraint(latitude=lambda cell: cell < maxlat)
    box_constraint = (longitude_constraint1 & longitude_constraint2 &
                      latitude_constraint1 & latitude_constraint2)
    return box_constraint
