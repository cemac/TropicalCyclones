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
import matplotlib.pyplot as plt  # KEEP ME HERE!!!


def checker(dd, hr, d0, dN, t0, tN):
    """checker
    Description:
        Check in correct time critera depending on set start and end time.
    Args:
        dd(int): day
        hr(int): hour
        d0(int): init day
        dN(int): end day
        t0(int): init hour
        tN(int): end hour
    Return:
        lead_time (int) or lead_time .false.
    """
    if dd < d0 or (dd == d0 and hr < t0 + 1):
        lead_time = False
    if dd > dN or (dd > dN - 1 and hr > tN - 1):
        lead_time = False
    lead_time = (dd - d0) * 24 + (hr - t0)
    return lead_time


def find_subplot_dims(n):
    """find_subplot_dims
    Description:
        Automatic format subplots for given number of plots
    Args:
        n (int): number of ensemble members
    Return:
        n_rows (int): number of rows
        n_cols (int): number of columns
    """
    n_cols = np.ceil(n**0.5).astype(int)
    n_rows = np.ceil(1. * n / n_cols).astype(int)
    return n_rows, n_cols


def uv(df, u_cons, t_cons, p_cons, b_cons):
    """uv
    Description:
        Extract velocity field under constants
    Args:
        df (str): file name
        u_cons (iris constraint): velocity constraint
        t_cons (iris constraint): temporal constraint (set by time_const)
        p_cons (iris constraint): pressure constraint
        b_cons (iris constraint): box constraint (set by box_constraint)
    Return:
        u (iris cube): velocity field matching constraints
    """
    with iris.FUTURE.context(cell_datetime_objects=True):
        u = iris.load_cube(df, u_cons).extract(t_cons & p_cons & b_cons)
    return u


def winds(u, v, em):
    """winds
    Description:
        Create cube of WindSpeeds for all ensemble members
    Args:
        u (iris cube): zonal velocity cube (produced by uv)
        v (iris cube): meridional velocity cube (produced by uv)
        em (int): ensemble member number
    Return:
        u (iris cube): zonal velocity cube with assigned em member
        v (iris cube): meridional velocity cube with assigned em member
        ws (iris cube): windspeed cube with assigned em member
    """
    ws = (u**2 + v**2)**0.5
    coord = iris.coords.AuxCoord(em, long_name='ensemble_member')
    u.add_aux_coord(coord)
    v.add_aux_coord(coord)
    ws.add_aux_coord(coord)
    return u, v, ws


def annotate(axs, StringFormat, xy):
    """annotate
    Description:
        Create cube of WindSpeeds for all ensemble members
    Args:
        axs (fig axes): figure axes
        StringFormat (str): Regex string
        xy:
    Return:
        Adds annoation to axs
    """
    # Add initial time, valid time etc.
    axs.annotate(StringFormat, xy=xy, xycoords='figure fraction',
                 horizontalalignment='right', verticalalignment='top',
                 color='k', fontsize=15)


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


def plot_wspeed(ax, ws):
    """plot_wspeed
    Description:
        Plot WindSpeed
    Args:
        axs (fig axes): figure axes
        ws (iris cube): windspeed cube
    Return:
        wspeed (figure): WindSpeed contourf plot
    """
    levels = np.arange(11) * 5
    cmap = mpl_cm.get_cmap('plasma')
    x = ws.coord('longitude').points
    y = ws.coord('latitude').points
    X, Y = np.meshgrid(x, y)
    wspeed = ax.contourf(X, Y, ws.data, cmap=cmap,
                         levels=levels, extend='both')
    return wspeed


def plot_winds(ax, u, v, mod):
    """plot_winds
    Description:
        Plot quiver of wind field
    Args:
        axs (fig axes): figure axes
        u (iris cube): windspeed cube
        v (iris cube): windspeed cube
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
    """map_formatter
    Description:
        Adds gridlines, countries and lat/lon labels to the plot.
    Args:
        axs (fig axes): figure axes
        params: preset but adjustable
    Return:
        fig with gridlines added
    """
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
    """find_centre_3h
    Description:
    Loads the track data for a specific time. This function is tailored to johns
    data, will have to be rewritten depending on how the users track data
    is stored. If you want just a single domain simply create a function that
    returns the centre of this domain.
    Args:
        md(int): day
        TT(int): hour
        em(int): ensemble number
        v_day(list): list of days
        v_time (list): list of times
        mod (str): model
        froot (str): file path
        fpat (str): regex filename pattern
        dom_r (float):
    Return:
        minlat, maxlat, minlon, maxlon: box domain
    """
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
    ir = fload.extract(iris.Constraint(longitude=lambda cell: minlon < cell <
                                       maxlon, latitude=lambda cell: minlat < cell < maxlat))
    return ir


def load_ens_members(em, fpath, x0, y0, data_constraint1, data_constraint2,
                     p_constraint):
    """load_ens_members
    Description:
    Args:
        em (int): ensemble number
        fpath (str): file path
        x0 (1D array):
        y0 (1D array):
        data_constraint1 (iris constraint): u_constraint
        data_constraint2 (iris constraint): v_constraint
        p_constraint (iris constraint): pressure constraint
    Return:
        v_azi_all: tangential azimuth velocity
        u_rad_all: radial azimuth velocity
        vert: vertical velocity
    """
    phi_interval = np.pi / 8
    ranges = np.arange(0, 500, 5)
    phis = np.arange(0, 2 * np.pi, phi_interval)
    df = fpath + '{0:02d}.pp'.format(em)
    u = iris.load_cube(df, data_constraint1).extract(p_constraint)
    v = iris.load_cube(df, data_constraint2).extract(p_constraint)

    i = 0

    for u_slc, v_slc in zip(u.slices(['latitude', 'longitude']), v.slices(['latitude', 'longitude'])):
        if i < 2:
            minlat = y0[0] - 5
            maxlat = y0[0] + 5
            minlon = x0[0] - 5
            maxlon = x0[0] + 5
        else:
            minlat = y0[i / 2 - 1] - 5
            maxlat = y0[i / 2 - 1] + 5
            minlon = x0[i / 2 - 1] - 5
            maxlon = x0[i / 2 - 1] + 5

        b_constraint = box_constraint(minlat, maxlat, minlon, maxlon)
        u_box = u_slc.extract(b_constraint)
        v_box = u_slc.extract(b_constraint)
        vrt = calc_vrt_spherical(u_box, v_box)
        cenlat, cenlon, max_val = max_vals(vrt)
        for r in ranges:
            for phi in phis:
                xpoi = cenlon + 0.009 * r * np.cos(phi)  # Check the 0.009
                ypoi = cenlat + 0.009 * r * np.sin(phi)
                new_point = [('latitude', ypoi), ('longitude', xpoi)]
                newu = u_slc.interpolate(
                    new_point, iris.analysis.Linear()).data
                newv = v_slc.interpolate(
                    new_point, iris.analysis.Linear()).data
                aziv = -newu * np.sin(phi) + newv * np.cos(phi)
                radu = newu * np.cos(phi) + newv * np.sin(phi)
                if phi == 0:
                    vazi = aziv
                    urad = radu
                else:
                    vazi = np.append(vazi, aziv)
                    urad = np.append(urad, radu)

            if r == 0:
                Vazi = np.mean(vazi)
                Urad = np.mean(urad)
            else:
                Vazi = np.append(Vazi, np.mean(vazi))
                Urad = np.append(Urad, np.mean(urad))
        if i == 0:
            v_azi_all = Vazi
            u_rad_all = Urad
        else:
            v_azi_all = np.dstack((v_azi_all, Vazi))
            u_rad_all = np.dstack((u_rad_all, Urad))

        i = i + 1

    return v_azi_all, u_rad_all


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


def calc_vrt_spherical(u, v):
    """calc_vrt_spherical
    Description:
        Calculates vorticity in spherical coordinates using the following:
        omega = (1/(a cos(lat)))*(dv/dlon - d(u cos(lat))/dlat)
    Args:
        u (iris cube):
        v (iris cube):
    Returns:
        vrt: vorticity
    """
    lats = u.coord('latitude').points
    lons = u.coord('longitude').points
    LONS, LATS = np.meshgrid(lons, lats)

    a = 6371000.  # Radius of Earth
    ddeg = abs(lats[1] - lats[0])
    drad = np.pi * ddeg / 180.  # Grid spacing (assumed constant)

    # Calculate dv/dlon, i.e. for each latitude
    for i in np.arange(len(u.data)):
        dvdlon = calc_grad(v.data[i, :], drad)
        if i == 0:
            v_lon = dvdlon
        else:
            v_lon = np.vstack((v_lon, dvdlon))

    lats_rad = LATS * np.pi / 180.
    ucoslat = u.data * np.cos(lats_rad)

    # Calculate d(u sin(lat))/dlat, i.e. for each lon
    for i in np.arange(len(u.data[0])):
        dudlat = calc_grad(ucoslat[:, i], drad)
        if i == 0:
            u_lat = dudlat
        else:
            u_lat = np.vstack((u_lat, dudlat))
    u_lat = np.swapaxes(u_lat, 0, 1)

    vrt = u[:]
    vrt.data = (1 / (a * np.cos(lats_rad))) * (v_lon - u_lat)
    vrt.units = 's-1'
    vrt.standard_name = 'atmosphere_relative_vorticity'
    vrt.long_name = 'Calculated using spherical coordinates finite difference'
    vrt.attributes = {
        'source': 'Calculated using (u,v) from Met Office Unified Model', 'um_version': '10.6'}
    return vrt


def plot_hovmoller(v_azi, vrad, outfile, em):
    """plot_hovmoller
    Description
    Args:
        v_azi: tangential azimuth velocity
        vrad: radial azimuth velocity
        outfile (str):
        em (str):
    Return:
        hovmoller plot
    """
    ranges = np.arange(0, 500, 5)
    fig = plt.figure(figsize=(16, 16))
    data = v_azi[0]
    data = np.swapaxes(data, 0, 1)
    times = np.arange(41)
    fig = plt.figure(1)
    ax = fig.add_subplot(1, 3, 1)
    ax.set_xlabel('Radius (km)', fontsize=18)
    ax.set_ylabel('Forecast time', fontsize=18)
    hovmol = ax.contourf(ranges, times, data, cmap='viridis', extend='both')
    # Contour mean tangential wind
    cbar = plt.colorbar(hovmol, orientation='horizontal', extend='both',
                        fraction=0.046, pad=0.09)
    cbar.set_label('dAzimuthal velocity (ms$^{-1}$)', size=14)
    cbar.ax.tick_params(labelsize=14)
    ax = fig.add_subplot(1, 3, 2)
    data = vrad[0]
    data = np.swapaxes(data, 0, 1)
    ax.set_xlabel('Radius (km)', fontsize=18)
    ax.set_ylabel('Forecast time', fontsize=18)
    hovmol = ax.contourf(ranges, times, data, cmap='viridis', extend='both')
    # Contour mean tangential wind
    cbar = plt.colorbar(hovmol, orientation='horizontal', extend='both',
                        fraction=0.046, pad=0.09)
    cbar.set_label('Radial velocity (ms$^{-1}$)', size=14)
    cbar.ax.tick_params(labelsize=14)
    fig.suptitle('Simulation em' + str(em))
    plt.tightlayout()
    plt.savefig(outfile)
    plt.close()
