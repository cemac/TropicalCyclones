# -*- coding: utf-8 -*-
"""ToolKit

.. module:: ToolKit
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
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.ticker as mticker
import matplotlib.cm as mpl_cm


def checker(dd, hr, d0, dN, t0, tN):
    if dd < d0 or (dd == d0 and hr < t0 + 1):
        lead_time = False
    if dd > dN or (dd > dN - 1 and hr > tN - 1):
        lead_time = False
    lead_time = (dd - d0) * 24 + (hr - t0)
    print('Plotting data for dd = {0}, hr = {1}'.format(dd, hr))
    return lead_time


def find_subplot_dims(n):
    n_cols = np.ceil(n**0.5).astype(int)
    n_rows = np.ceil(1. * n / n_cols).astype(int)
    return n_rows, n_cols


def uv(df, u_cons, t_cons, p_cons, b_cons):
    with iris.FUTURE.context(cell_datetime_objects=True):
        u = iris.load_cube(df, u_cons).extract(t_cons & p_cons & b_cons)
    return u


def winds(u, v, em):
    ws = (u**2 + v**2)**0.5
    coord = iris.coords.AuxCoord(em, long_name='ensemble_member')
    u.add_aux_coord(coord)
    v.add_aux_coord(coord)
    ws.add_aux_coord(coord)
    return u, v, ws


def annotate(axs, StringFormat, xy):
    # Add initial time, valid time etc.
    axs.annotate(StringFormat, xy=xy, xycoords='figure fraction',
                 horizontalalignment='right', verticalalignment='top',
                 color='k', fontsize=15)


def label_fixer(i, ncols, nrows):
    # Make sure only the outside axis are labelled
    if i % ncols != 0:
        left_label = False
    else:
        left_label = True
    if i < ncols * (nrows - 1):
        bottom_label = False
    else:
        bottom_label = True
    return left_label, bottom_label


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
                  left_label=True, res='10m'):
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


def find_centre_3h(md, TT, em, v_day, v_time, mod, froot, fpat, dom_r):
    '''
    Loads the track data for a specific time. This function is tailored to my
    own data, will have to be rewritten depending on how the users track data
    is stored. If you want just a single domain simply create a function that
    returns the centre of this domain.
    '''
    if mod == 'glm':
        fpatfull = froot + '_newinterp_' + fpat
    else:
        fpatfull = froot + fpat
    filepath = fpatfull.format(md, TT, em, mod)
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
    try:
        minlat = cenlat[0] - dom_r
        maxlat = cenlat[0] + dom_r
        minlon = cenlon[0] - dom_r
        maxlon = cenlon[0] + dom_r
    except IndexError:
        return False
    return minlat, maxlat, minlon, maxlon


def box_constraint(minlat, maxlat, minlon, maxlon):
    # Create constraint to extract data from cube over a certain region
    longitude_constraint1 = iris.Constraint(
        longitude=lambda cell: cell > minlon)
    longitude_constraint2 = iris.Constraint(
        longitude=lambda cell: cell < maxlon)
    latitude_constraint1 = iris.Constraint(latitude=lambda cell: cell > minlat)
    latitude_constraint2 = iris.Constraint(latitude=lambda cell: cell < maxlat)
    tc_box_constraint = (longitude_constraint1 & longitude_constraint2
                         & latitude_constraint1 & latitude_constraint2)
    return tc_box_constraint


def calculate_pot_temp(pressure, temperature):
    return temperature * (1000 / pressure)**(0.286)


def calc_grad(data, dx):
    # Basic central differencing
    # Data should be 1-d
    grad = np.zeros(data.shape)
    grad[0] = (1 / dx) * (data[1] - data[0])
    grad[-1] = (1 / dx) * (data[-1] - data[-2])
    grad[1:-1] = (1 / (2 * dx)) * (data[2:] - data[:-2])
    return grad


def extracter(fload, minlon, maxlon, minlat, maxlat):
    ir = fload.extract(iris.Constraint(longitude=lambda cell: minlon < cell
                       < maxlon, latitude=lambda cell: minlat < cell < maxlat))
    return ir
