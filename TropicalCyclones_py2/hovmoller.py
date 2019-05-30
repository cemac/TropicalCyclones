# -*- coding: utf-8 -*-
"""Diags

.. module:: Hovmoller Plotter
    :platform: Unix
    :synopis: Plot Hovmoller of azimuthal velocity

.. moduleauthor: Sam Hardy(UoL), John Ashcroft (UoL), Helen Burns (CEMAC-UoL)

.. description: This module was developed by CEMAC as part of the WCSSP
                Project.

   :copyright: © 2019 University of Leeds.
   :license: MIT

Example:
    To use::

Memebers:

.. CEMAC_TropicalCyclones:
   https://github.com/cemac/TropicalCyclones
"""
from __future__ import print_function
import sys
import iris
import iris.analysis.calculus
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import toolkit as tct
mpl.pyplot.switch_backend('Agg')


def set_ofile():
    """set_ofile
    Description:
        Sets outputfile locationa and name
    Args: None
    Returns:
        ofile (str): path for output
    """
    ofile = ('plots/hovs/{0:04d}_{1:02d}Zhovmoller_4p4_{0:04d}_{1:02d}Z_'
             + 'em{2:02d}.png')
    return ofile


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


def calc_azimuth_vels(ens, fpath, x_0, y_0, constraints):
    """calc_azimuth_vels
    Description:
     Calculate a azimuthal velocitiesfor a give ensemble member, give data
     location and storm box. Will save the numpy arrays in order to save time
     in the future.
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
        vort: vorticty
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

        b_constraint = tct.box_constraint(minlat, maxlat, minlon, maxlon)
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
                    vrtcyl = newvrt
                else:
                    vazi = np.append(vazi, aziv)
                    urad = np.append(urad, radu)
                    vrtcyl = np.append(vrtcyl, newvrt)
            if num == 0:
                vrt_cyl = np.mean(vrtcyl)
                v_azi = np.mean(vazi)
                u_rad = np.mean(urad)
            else:
                vrt_cyl = np.append(vrt_cyl, np.mean(vrtcyl))
                v_azi = np.append(v_azi, np.mean(vazi))
                u_rad = np.append(u_rad, np.mean(urad))
        if i == 0:
            v_azi_all = v_azi
            u_rad_all = u_rad
            vrt_all = vrt_cyl
        else:
            v_azi_all = np.dstack((v_azi_all, v_azi))
            u_rad_all = np.dstack((u_rad_all, u_rad))
            vrt_all = np.dstack((vrt_all, vrt_cyl))

        i = i + 1
    np.save('v_azi_EM{0:02d}.npy'.format(ens), v_azi_all)
    np.save('u_rad_EM{0:02d}.npy'.format(ens), u_rad_all)
    np.save('vrt_all_EM{0:02d}.npy'.format(ens), vrt_all)
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
        vrt (array): vorticity
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
        v_azi(T,Y) (array): tangential azimuth velocity
        vrad(T,Y) (array): radial azimuth velocity
        vrt(T,Y) (array): Vorticity
        outfile (str): path and file name for plot
        ens (str): ensemble member number
    Return:
        hovmoller plot
    """
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    ranges = np.arange(0, 5 * (v_azi[0]).shape[0], 5)
    times = np.arange((v_azi[0]).shape[1])
    fig = plt.figure(figsize=(16, 16))
    plot_vtan_rad(v_azi[0], v_azi[0], ['dAzimuthal velocity', '(ms$^{-1}$)'],
                  fig, [ranges, times], 1)
    plot_vtan_rad(vrad[0], v_azi[0], ['Radial velocity', '(ms$^{-1}$)'],
                  fig, [ranges, times], 2)
    plot_vort(vrt[0], fig, ranges, times, 3)
    fig.suptitle('Simulation EM{0:02d}'.format(ens), fontsize=20)
    plt.tight_layout()
    fig.subplots_adjust(top=0.9)
    plt.savefig(outfile)
    plt.close()


def plot_vtan_rad(data, v_azi, var_units, fig, y_t, i, cmap='PuOr_r',
                  tend='Y'):
    """plot_vtan_rad
    Description:
        Plots a subplot i on fig axis fig of data with contours of v_azi in
        Y-T space defaulting to plotting tendency
    Args:
        data(T, Y) (2D array): data to be plotted e.g radial azimuth velocity
        Vazi(T, Y) (2D array): idealy tangential azimuth velocity
        var_units (list):
            var (str): Variable name
            units (str): units
        fig (mpl axis): axes to plot to
        y_t (list): list of dimensions
            Y (1D array): range (radius (km))
            T (1D array): time (hours)
        i (int): suplot number
        cmap (mpl colourmap): Choose colour map default = PuOr_r
        tend (str): Plot tendency (Y) or not
    Return:
        hovmoller plot
    """
    fig_letter = ['a) ', 'b) ', 'c) ']
    data = np.swapaxes(data, 0, 1)
    v_azi = np.swapaxes(v_azi, 0, 1)
    axs = fig.add_subplot(1, 3, i)
    axs.set_xlabel('Radius (km)', fontsize=18)
    axs.set_ylabel('Forecast time (h)', fontsize=18)
    if tend == 'Y':
        d_tend = (data[1::, :] - data[0:-1, :]) / 2
        clab = 'Tendency '
    else:
        d_tend = data[1::]
        clab = ''
    maxi = int(np.array(abs(d_tend.max()), abs(d_tend.min())).max())
    step = maxi / 8.0
    fill = axs.contourf(y_t[0], y_t[1][1::], d_tend,
                        levels=np.arange(-maxi, maxi, step) + step / 2.0,
                        cmap=cmap, extend='both')
    axs.contour(y_t[0], y_t[1][1::], v_azi[1::, :],
                levels=[np.arange(int(v_azi.min()), int(v_azi.max()), 10)],
                colors='black', linewidths=3)
    axs.contour(y_t[0], y_t[1][1::], d_tend, levels=[0], colors='grey',
                linewidths=2)
    # Contour mean tangential wind
    tct.annotate(axs, fig_letter[i] + var_units[0] + ' ' + var_units[1],
             [i * 0.3, 0.85])
    cbar = plt.colorbar(fill, orientation='horizontal', extend='both',
                        fraction=0.046, pad=0.09)
    cbar.set_ticks([np.arange(-maxi, maxi, int(step * 2.0))])
    cbar.set_label(clab + var_units[0], size=18)
    cbar.ax.tick_params(labelsize=18)


def plot_vort(data, fig, ranges, times, i, cmap='PuOr_r'):
    """plot_vtan_rad
    Description:
        Plots a subplot vorticity i on fig axis fig of data with contours of
        v_azi in Y-T space defaulting to plotting tendency
    Args:
        data(T, Y) (2D array): data to be plotted e.g radial azimuth velocity
        fig (mpl axis): axes to plot to
        y_t (list): list of dimensions
        ranges (1D array): range (radius (km))
        times (1D array): time (hours)
        i (int): suplot number
        cmap (mpl colourmap): Choose colour map default = PuOr_r
    Return:
        hovmoller plot
    """
    axs = fig.add_subplot(1, 3, i)
    data = np.swapaxes(data, 0, 1)
    d_tend = 1000 * (data[1::, :] - data[0:-1, :]) / 2
    axs.set_xlabel('Radius (km)', fontsize=18)
    axs.set_ylabel('Forecast time (h)', fontsize=18)
    maxi = int(np.array(abs(d_tend.max()), abs(d_tend.min())).max())
    step = maxi / 8.0
    fill = axs.contourf(ranges, times[1::], d_tend, cmap=cmap,
                        levels=np.arange(-maxi, maxi, step) + step / 2.0,
                        extend='both')
    axs.contour(ranges, times[1::], d_tend, levels=[0], colors='grey',
                linewidths=2)
    tct.annotate(axs, 'c) Vorticity', [i * 0.3, 0.85])
    # Contour Vorticity
    cbar = plt.colorbar(fill, orientation='horizontal', extend='both',
                        fraction=0.046, pad=0.09)
    cbar.set_label('Tendency Vorticity x $10 ^{-3}$', size=18)
    cbar.ax.tick_params(labelsize=18)
    cbar.set_ticks([np.arange(-maxi, maxi, step * 2.0)])
    cbar.ax.set_yticklabels(np.arange(int(data.min()), int(data.max()), 0.002),
                            fontsize=16, weight='bold')


class HovPlotter(object):
    '''Hovmoller Plotter
        Plots Hovmoller using configuration file

    Members:
        hovloop: loop over ensemble members at give time
        plot_ens: plot WindSpeed plot for each ensemble member
        hovplotter: call hovoller plotter for time and memeber
    '''

    def __init__(self, configfile='../configfiles/configfile.csv',
                 stashfile='../configfiles/stashvars'):
        '''
        Args:
            configfile (string): filepath to configuration settings
            stashfile (string): filepath to shashfile codes
        '''
        # Read configuration file to import settings
        conf_df = pd.read_csv(configfile + '.csv')
        conf_df = conf_df.set_index('var:').transpose()
        # Define all constraints and locate data.
        self.data_loc = (conf_df.data_root[0] + '/' + conf_df.Storm[0] + '/' +
                         conf_df.suitno[0] + '/' + conf_df.run[0] + '/' +
                         conf_df.mmdd_hr[0] + '/' + conf_df.vars[0])
        self.data_name = conf_df.data_name[0]
        self.ens_members = np.arange(int(conf_df.ens_members[0]))
        self.plev = int(conf_df.plev[0])
        # dates
        self.dates = [int(conf_df.yr[0]), int(conf_df.mth[0]),
                      int(conf_df.md[0]), int(conf_df.TT[0])]
        self.track = [float(conf_df.domain_rad[0]), conf_df.track_data_root[0],
                      conf_df.track_data_metadata[0]]
        print('Loaded configuration settings')
        # add stash codes
        self.stashfile = stashfile

    def stash_vars(self):
        """stash_vars
        Description:
            Uses Stash file to create iris cube contraints (variables)
        Args: none
        Returns:
            u_constraint (iris cube constraint): Uvel constraint
            v_constraint (iris cube constraint): Vvel constraint
            p_constraint (iris cube constraint): Pressure constraint
        """
        stash_df = pd.read_csv(self.stashfile + '.csv')
        u_stash = stash_df.u_stash[0]
        v_stash = stash_df.v_stash[0]
        u_constraint = iris.AttributeConstraint(STASH=u_stash)
        v_constraint = iris.AttributeConstraint(STASH=v_stash)
        p_constraint = iris.Constraint(pressure=self.plev)
        return u_constraint, v_constraint, p_constraint

    def hovloop(self):
        """hovloop
        Description:
            Single plots for each ensemble member of tangential, radial and
            vertical wind speed. Allows to loop over a number of dates and
            ensemble members
        Args:
            none
        Returns:
            Hovmoller plot for a day and time for a set ensemble member
        """
        mmdd_time = self.dates
        mmdd, time = mmdd_time[-2], mmdd_time[-1]
        for mmdd in [mmdd]:
            for hour in [time]:
                for ens in self.ens_members:
                    self.hovplotter(mmdd, hour, ens)

    def hovplotter(self, mmdd, hour, ens):
        """hovplotter
        Description:
        Args:
            mmdd (int): month-day MMDD
            hour (int): hour
            ens (int): ensemble member
        Returns:
            hovmoller plot of ensemble memeber for set time
        """
        ofile = set_ofile()
        fpath = self.data_loc.format(
            mmdd, hour) + self.data_name.format(mmdd, hour)
        outfile = ofile.format(mmdd, hour, ens)
        constraints = self.stash_vars()
        y_0, x_0 = self.storm_center(mmdd, hour, ens)
        vtan, vrad, vrt = calc_azimuth_vels(ens, fpath, x_0, y_0,
                                            constraints)
        plot_hovmoller(vtan, vrad, vrt, outfile, ens)

    def storm_center(self, mmdd, hour, ens):
        """storm_center
        Description:
            Extracts storm center information from stored track data
        Args:
            mmdd (int): month-day MMDD
            hour (int): hour
            ens (int): ensemble member
        Returns:
            x_0 centre longitude
            y_0 center latitue
        """
        filepath = str(self.track[1]) + \
            '_4p4_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(mmdd, hour, ens)
        try:
            track_data = np.load(filepath)
        except ValueError:
            print('Expected npz format')
            print('Need storm center info')
            sys.exit(0)
        lats = track_data['lats']
        lons = track_data['lons']
        [y_0, x_0] = [lats, lons]  # centre of storm
        return y_0, x_0
