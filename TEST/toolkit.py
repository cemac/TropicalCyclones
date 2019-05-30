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
from __future__ import print_function
import iris
import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib
import matplotlib.ticker as mticker
import matplotlib.cm as mpl_cm
import matplotlib.pyplot as plt  # KEEP ME HERE!!!
matplotlib.pyplot.switch_backend('Agg')


def checker(day, hour, times):
    """checker
    Description:
        Check in correct time critera depending on set start and end time.
    Args:
        day(int): day
        hour(int): hour
        times (list): list of start and end times
            day_0(int): init day
            day_n(int): end day
            time_0(int): init hour
            time_n(int): end hour
    Return:
        lead_time (int) or lead_time .false.
    """
    day_0, day_n, time_0, time_n = times
    if day < day_0 or (day == day_0 and hour < time_0 + 1):
        lead_time = False
    if day > day_n or (day > day_n - 1 and hour > time_n - 1):
        lead_time = False
    lead_time = (day - day_0) * 24 + (hour - time_0)
    return lead_time


def find_subplot_dims(n_ensembles):
    """find_subplot_dims
    Description:
        Automatic format subplots for given number of plots
    Args:
        n (int): number of ensemble members
    Return:
        n_rows (int): number of rows
        n_cols (int): number of columns
    """
    n_cols = np.ceil(n_ensembles**0.5).astype(int)
    n_rows = np.ceil(1. * n_ensembles / n_cols).astype(int)
    return n_rows, n_cols


def vel_extract(d_file, u_cons, t_cons, p_cons, b_cons):
    """uv
    Description:
        Extract velocity field under constants
    Args:
        d_file (str): file name
        u_cons (iris constraint): velocity constraint
        t_cons (iris constraint): temporal constraint (set by time_const)
        p_cons (iris constraint): pressure constraint
        b_cons (iris constraint): box constraint (set by box_constraint)
    Return:
        u (iris cube): velocity field matching constraints
    """
    with iris.FUTURE.context(cell_datetime_objects=True):
        vel = iris.load_cube(d_file, u_cons).extract(t_cons & p_cons & b_cons)
    return vel


def winds(uvel, vvel, ens):
    """winds
    Description:
        Create cube of WindSpeeds for all ensemble members
    Args:
        uvel (iris cube): zonal velocity cube (produced by uv)
        vvel (iris cube): meridional velocity cube (produced by uv)
        ens (int): ensemble member number
    Return:
        uvel (iris cube): zonal velocity cube with assigned em member
        vvel (iris cube): meridional velocity cube with assigned em member
        w_speed (iris cube): windspeed cube with assigned em member
    """
    w_speed = (uvel**2 + vvel**2)**0.5
    coord = iris.coords.AuxCoord(ens, long_name='ensemble_member')
    uvel.add_aux_coord(coord)
    vvel.add_aux_coord(coord)
    w_speed.add_aux_coord(coord)
    return uvel, vvel, w_speed


def annotate(axs, str_format, xy):
    """annotate
    Description:
        Create cube of WindSpeeds for all ensemble members
    Args:
        axs (fig axes): figure axes
        str_format (str): Regex string
        xy:
    Return:
        Adds annoation to axs
    """
    # Add initial time, valid time etc.
    bbox_args = dict(boxstyle="round", fc="0.8")
    axs.annotate(str_format, xy=xy,
                 xycoords='figure fraction', xytext=(40, 20),
                 textcoords='offset points', ha="right", va="top",
                 bbox=bbox_args, fontsize=16)


def label_fixer(i, ncols, nrows):
    """label_fixer
    Description:
        Create cube of WindSpeeds for all ensemble members
    Args:
        i (int): plot number (iterable)
        ncols (int): number of columns from find_subplot_dims
        nrows (int): number of rows from find_subplot_dims
    Return:
        left_label (str): left labels
        bottom_label (str): bottom labes
    """
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


def plot_wspeed(ax_ws, w_speed):
    """plot_wspeed
    Description:
        Plot WindSpeed
    Args:
        ax_ws (fig axes): figure axes
        w_speed (iris cube): windspeed cube
    Return:
        wspeed (figure): WindSpeed contourf plot
    """
    levels = np.arange(11) * 5
    cmap = mpl_cm.get_cmap('plasma')
    x = w_speed.coord('longitude').points
    y = w_speed.coord('latitude').points
    x, y = np.meshgrid(x, y)
    wspeed = ax_ws.contourf(x, y, w_speed.data, cmap=cmap,
                            levels=levels, extend='both')
    return wspeed


def plot_winds(ax_ws, uvel, vvel, mod):
    """plot_winds
    Description:
        Plot quiver of wind field
    Args:
        ax_ws (fig axes): figure axes
        uvel (iris cube): windspeed cube
        vvel (iris cube): windspeed cube
        mod (str): model
    Return:
        wspeed (figure)
    """
    if mod == '4p4':
        scale = 300.
        step = 12
    else:
        scale = 400.
        step = 8  # Need to trial and error these.

    x = uvel.coord('longitude').points
    y = uvel.coord('latitude').points
    x, y = np.meshgrid(x, y)

    arrows = ax_ws.quiver(x[::step, ::step], y[::step, ::step],
                          uvel.data[::step, ::step],
                          vvel.data[::step, ::step], angles='xy', scale=scale,
                          transform=ccrs.PlateCarree())
    if uvel.coord('ensemble_member').points[0] == 0:
        ax_ws.quiverkey(arrows, 0.85, 0.1, 20, '20ms$^{-1}$',
                        coordinates='figure', fontproperties={'size': 15})


def map_formatter(var_ax, tick_base_x=15.0, tick_base_y=15.0, labelsize=20,
                  top_label=False, bottom_label=True, right_label=False,
                  left_label=True, res='10m'):
    """map_formatter
    Description:
        Adds gridlines, countries and lat/lon labels to the plot.
    Args:
        var_ax (fig axes): figure axes
        params: preset but adjustable
    Return:
        fig with gridlines added
    """
    grl = var_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                           linewidth=0.75, color='k', linestyle=':')
    grl.xlabels_top = top_label
    grl.ylabels_right = right_label
    grl.ylabels_left = left_label
    grl.xlabels_bottom = bottom_label
    grl.xlocator = mticker.MultipleLocator(base=tick_base_x)
    grl.ylocator = mticker.MultipleLocator(base=tick_base_y)
    var_ax.coastlines(resolution=res, color='k', linewidth=1)
    grl.xlabel_style = {'size': labelsize}
    grl.ylabel_style = {'size': labelsize}
    grl.xformatter = LONGITUDE_FORMATTER
    grl.yformatter = LATITUDE_FORMATTER


def find_centre_3h(info, finfo):
    """find_centre_3h
    Description:
    Loads the track data for a specific time. This function is tailored to
    johns data, will have to be rewritten depending on how the users track data
    is stored. If you want just a single domain simply create a function that
    returns the centre of this domain.
    Args:
        info (list): metadata
            mmdd(int): day
            time(int): hour
            em(int): ensemble number
            mod (str): model
        finfo (list): file info
            dom_r (float):
            froot (str): file path
            fpat (str): regex filename pattern
    Return:
        minlat, maxlat, minlon, maxlon: box domain
    """
    if info[-1] == 'glm':
        fpatfull = finfo[1] + '_newinterp_' + finfo[2]
    else:
        fpatfull = finfo[1] + finfo[2]
    filepath = fpatfull.format(info[0], info[2], info[3], info[4])
    track_data = np.load(filepath)
    lats = track_data['lats']
    lons = track_data['lons']
    days = track_data['vt_days']
    times = track_data['vt_times']
    index1 = np.where(days == info[1])
    lats = lats[index1]
    lons = lons[index1]
    times = times[index1]
    index2 = np.where(times == info[2])
    cenlat = lats[index2]
    cenlon = lons[index2]
    try:
        min_max = [cenlat[0] - finfo[0], cenlat[0] + finfo[0],
                   cenlon[0] - finfo[0], cenlon[0] + finfo[0]]
    except IndexError:
        return False
    return min_max


def box_constraint(minlat, maxlat, minlon, maxlon):
    """box_constraint
    Description:
        Contrain cube using lat lons
    Args:
        minlat (int): minium latitue
        minlon (int): minimum longitude
        maxlat (int): maximum latitue
        maxlon (int): maximum longitude
    Return:
        tc_box_constraint (iris Constraints)
    """
    # Create constraint to extract data from cube over a certain region
    longitude_constraint1 = iris.Constraint(
        longitude=lambda cell: cell > minlon)
    longitude_constraint2 = iris.Constraint(
        longitude=lambda cell: cell < maxlon)
    latitude_constraint1 = iris.Constraint(latitude=lambda cell: cell > minlat)
    latitude_constraint2 = iris.Constraint(latitude=lambda cell: cell < maxlat)
    tc_box_constraint = (longitude_constraint1 & longitude_constraint2 &
                         latitude_constraint1 & latitude_constraint2)
    return tc_box_constraint


def calculate_pot_temp(pressure, temperature):
    """calculate_pot_temp
    Description:
        Temp * (1000/Pressure)^0.286
    Return potentail temp from pressure and temperature
    """
    return temperature * (1000 / pressure)**(0.286)


def calc_grad(data, dx):
    """calc_grad
    Description:
        Basic central differencing
    Args:
        data (1D array): 1D data
        dx: spacing
    Return:
        gradient (1D array)
    """
    # Basic central differencing
    # Data should be 1-d
    grad = np.zeros(data.shape)
    grad[0] = (1 / dx) * (data[1] - data[0])
    grad[-1] = (1 / dx) * (data[-1] - data[-2])
    grad[1:-1] = (1 / (2 * dx)) * (data[2:] - data[:-2])
    return grad


def extracter(fload, minlon, maxlon, minlat, maxlat):
    """extracter
    Description:
    Args:
        fload (iris cube): loaded file
        minlon (int): minimum longitude
        maxlon (int): maximum longitude
        minlat (int): minimum latitude
        maxlat (int): maximum latitude
    Return:
        ir (iris cube): Contrained iris cube
    """
    cube = fload.extract(iris.Constraint(longitude=lambda cell: minlon < cell <
                                         maxlon, latitude=lambda cell: minlat
                                         < cell < maxlat))
    return cube


def load_ens_members(ens, fpath, x_0, y_0, constraints):
    """load_ens_members
    Description:
    Args:
        ens (int): ensemble number
        fpath (str): file path
        x_0 (1D array):
        y_0 (1D array):
        contraints (list):
            u_constraint (iris constraint): u_constraint
            v_constraint (iris constraint): v_constraint
            p_constraint (iris constraint): pressure constraint
    Return:
        v_azi_all: tangential azimuth velocity
        u_rad_all: radial azimuth velocity
        vert: vertical velocity
    """
    phi_interval = np.pi / 8
    ranges = np.arange(0, 500, 5)
    phis = np.arange(0, 2 * np.pi, phi_interval)
    d_file = fpath + '{0:02d}.pp'.format(ens)
    uvel = iris.load_cube(d_file, constraints[0]).extract(constraints[2])
    vvel = iris.load_cube(d_file, constraints[1]).extract(constraints[2])
    i = 0

    for u_slc, v_slc in zip(uvel.slices(['latitude', 'longitude']),
                            vvel.slices(['latitude', 'longitude'])):
        if i < 2:
            minlat = y_0[0] - 5
            maxlat = y_0[0] + 5
            minlon = x_0[0] - 5
            maxlon = x_0[0] + 5
        else:
            minlat = y_0[i / 2 - 1] - 5
            maxlat = y_0[i / 2 - 1] + 5
            minlon = x_0[i / 2 - 1] - 5
            maxlon = x_0[i / 2 - 1] + 5

        b_constraint = box_constraint(minlat, maxlat, minlon, maxlon)
        u_box = u_slc.extract(b_constraint)
        v_box = u_slc.extract(b_constraint)
        vrt = calc_vrt_spherical(u_box, v_box)
        cent_lls = max_vals(vrt)
        for num in ranges:
            for phi in phis:
                xpoi = cent_lls[1] + 0.009 * num * np.cos(phi)
                ypoi = cent_lls[0] + 0.009 * num * np.sin(phi)
                new_point = [('latitude', ypoi), ('longitude', xpoi)]
                newu = u_slc.interpolate(
                    new_point, iris.analysis.Linear()).data
                newv = v_slc.interpolate(
                    new_point, iris.analysis.Linear()).data
                aziv = -newu * np.sin(phi) + newv * np.cos(phi)
                radu = newu * np.cos(phi) + newv * np.sin(phi)
                newvrt = vrt.interpolate(
                    new_point, iris.analysis.Linear()).data
                if phi == 0:
                    vazi = aziv
                    urad = radu
                    vrtazi = newvrt
                else:
                    vazi = np.append(vazi, aziv)
                    urad = np.append(urad, radu)
                    vertazi = np.append(vrtazi, newvrt)
            if num == 0:
                vrt_azi = np.mean(vrtazi)
                v_azi = np.mean(vazi)
                u_rad = np.mean(urad)
            else:
                vrt_azi = np.append(vrt_azi, np.mean(vrtazi))
                v_azi = np.append(v_azi, np.mean(vazi))
                u_rad = np.append(u_rad, np.mean(urad))
        if i == 0:
            v_azi_all = v_azi
            u_rad_all = u_rad
            vrt_all = vrt_azi
        else:
            v_azi_all = np.dstack((v_azi_all, v_azi))
            u_rad_all = np.dstack((u_rad_all, u_rad))
            vrt_all = np.dstack((vrt_all, vrt_azi))

        i = i + 1
    return v_azi_all, u_rad_all, vrt_all


def max_vals(cube):
    """max_vals
    Description:
    Args:
        cube (iris cube): data
    Returns:
        cenlat: latitude coordinates of maximum value of data in cube.
        cenlon: longitude coordinates of maximum value of data in cube.
        maxval: maximum value of data in cube
    """
    # Find latitude/longitude coordinates of maximum value of data in cube.
    index = np.argmax(cube.data)
    indices = np.unravel_index(index, cube.data.shape)
    # Find the coordinates for this central position
    cenlat = cube.coord('latitude').points[indices[0]]
    cenlon = cube.coord('longitude').points[indices[1]]
    maxval = cube.data[indices]
    return cenlat, cenlon, maxval


def calc_vrt_spherical(uvel, vvel):
    """calc_vrt_spherical
    Description:
        Calculates vorticity in spherical coordinates using the following:
        omega = (1/(a cos(lat)))*(dv/dlon - d(u cos(lat))/dlat)
    Args:
        uvel (iris cube):
        vvel (iris cube):
    Returns:
        vrt: vorticity
    """
    lats = uvel.coord('latitude').points
    lons = uvel.coord('longitude').points
    ddeg = abs(lats[1] - lats[0])
    lons, lats = np.meshgrid(lons, lats)

    r_e = 6371000.  # Radius of Earth
    drad = np.pi * ddeg / 180.  # Grid spacing (assumed constant)
    # Calculate dv/dlon, i.e. for each latitude
    for i in np.arange(len(uvel.data)):
        dvdlon = calc_grad(vvel.data[i, :], drad)
        if i == 0:
            v_lon = dvdlon
        else:
            v_lon = np.vstack((v_lon, dvdlon))

    lats_rad = lats * np.pi / 180.
    ucoslat = uvel.data * np.cos(lats_rad)
    # Calculate d(u sin(lat))/dlat, i.e. for each lon
    for i in np.arange(len(uvel.data[0])):
        dudlat = calc_grad(ucoslat[:, i], drad)
        if i == 0:
            u_lat = dudlat
        else:
            u_lat = np.vstack((u_lat, dudlat))
    u_lat = np.swapaxes(u_lat, 0, 1)

    vrt = uvel[:]
    vrt.data = (1 / (r_e * np.cos(lats_rad))) * (v_lon - u_lat)
    vrt.units = 's-1'
    vrt.standard_name = 'atmosphere_relative_vorticity'
    vrt.long_name = 'Calculated using spherical coordinates finite difference'
    vrt.attributes = {'source':
                      'Calculated using (u,v) from Met Office Unified Model',
                      'um_version': '10.6'}
    return vrt


def plot_hovmoller(v_azi, vrad, vrt, outfile, ens):
    """plot_hovmoller
    Description
    Args:
        v_azi: tangential azimuth velocity
        vrad: radial azimuth velocity
        outfile (str):
        ens (str):
    Return:
        hovmoller plot
    """
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    ranges = np.arange(0, 500, 5)
    fig = plt.figure(figsize=(16, 16))
    data = v_azi[0]
    data = np.swapaxes(data, 0, 1)
    times = np.arange(41)
    fig = plt.figure(1)
    axs = fig.add_subplot(1, 3, 1)
    axs.set_xlabel('Radius (km)', fontsize=18)
    axs.set_ylabel('Forecast time (h)', fontsize=18)
    hovmol = axs.contourf(ranges, times, data, cmap='viridis', extend='both')
    # Contour mean tangential wind
    xy = [0.9, 0.1]
    annotate(axs, 'a) dAzimuthal velocity (ms$^{-1}$)', xy)
    cbar = plt.colorbar(hovmol, orientation='horizontal', extend='both',
                        fraction=0.046, pad=0.09)
    cbar.set_label('dAzimuthal velocity (ms$^{-1}$)', size=18)
    cbar.ax.tick_params(labelsize=18)
    axs = fig.add_subplot(1, 3, 2)
    data = vrad[0]
    data = np.swapaxes(data, 0, 1)
    axs.set_xlabel('Radius (km)', fontsize=18)
    axs.set_ylabel('Forecast time (h)', fontsize=18)
    hovmol = axs.contourf(ranges, times, data, cmap='viridis', extend='both')
    annotate(axs, 'b) Radial velocity (ms$^{-1}$)', xy)
    # Contour mean tangential wind
    cbar = plt.colorbar(hovmol, orientation='horizontal', extend='both',
                        fraction=0.046, pad=0.09)
    cbar.ax.tick_params(labelsize=18)
    axs = fig.add_subplot(1, 3, 3)
    data = vrt[0]
    data = np.swapaxes(data, 0, 1)
    times = np.arange(41)
    fig = plt.figure(1)
    axs.set_xlabel('Radius (km)', fontsize=18)
    axs.set_ylabel('Forecast time (h)', fontsize=18)
    hovmol = axs.contourf(ranges, times, data, cmap='viridis', extend='both')
    annotate(axs, 'c) Vorticity', xy)
    # Contour Vorticity
    cbar = plt.colorbar(hovmol, orientation='horizontal', extend='both',
                        fraction=0.046, pad=0.09)
    cbar.set_label('Vorticity', size=18)
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_yticklabels(np.arange(int(data.min()), int(data.max()), 0.002),
                            fontsize=16, weight='bold')
    fig = plt.gcf()
    fig.suptitle('Simulation EM' + str(ens), fontsize=20)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()
