# -*- coding: utf-8 -*-
"""Diags

.. module:: Diagnostics Plotter
    :platform: Unix
    :synopis: Wind Diagnostics

.. moduleauthor: John Ashcroft (UoL), Sam Hardy (UoL), Helen Burns (CEMAC-UoL)

..date: May 2019.

.. description: This module was developed by CEMAC as part of the WCSSP
                Project. Intial script improvements

   :copyright: © 2019 University of Leeds.
   :license: MIT.

Example:
    To use::

Memebers:

.. CEMAC_TropicalCyclones:
   https://github.com/cemac/TropicalCyclones
"""
from __future__ import print_function
import ast
import cPickle as pickle
import iris
import iris.analysis.calculus
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as image
import matplotlib.ticker as mticker
import toolkit as tct
mpl.pyplot.switch_backend('Agg')


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


class WS_Plotter(object):
    '''WS_Plotter
        Initialised extra varible and provides standard plots and
        plotting tools

    Members:
        windloop: loop over windspeed plots (ws_dayhour)
        ws_dayhour: loop over the day and hour and call stamp plotter
        plot_ens: plot WindSpeed plot for each ensemble member
        time_const: contrain for time and location
    '''

    def __init__(self, configfile='configfile2', stashfile='stashvars'):
        '''
        Args:
            configfile (string): filepath to configuration settings
            stashfile (string): filepath to shashfile codes
        '''
        # Read configuration file to import settings
        conf_df = pd.read_csv(configfile + '.csv')
        conf_df = conf_df.set_index('var:').transpose()
        # Define all constraints and locate data.
        self.data_loc = (conf_df.data_root[0] + '/' + conf_df.Storm[0] + '/'
                         + conf_df.suitno[0] + '/' + conf_df.run[0] + '/'
                         + conf_df.mmdd_hr[0] + '/' + conf_df.vars[0])
        self.data_name = conf_df.data_name[0]
        self.outfile_loc = conf_df.outfile_loc[0]
        self.outfile_name = conf_df.outfile_name.values
        self.ens_members = np.arange(int(conf_df.ens_members[0]))
        self.plev = int(conf_df.plev[0])
        self.model = conf_df.model[0]
        # Times
        init_day = int(conf_df.init_day[0])
        init_time = int(conf_df.init_time[0])
        final_day = int(conf_df.final_day[0])
        final_time = int(conf_df.final_time[0])
        self.times = [init_day, init_time, final_day, final_time]
        v_times = np.array(ast.literal_eval(conf_df.v_times[0]))
        v_days = np.arange(init_day, init_day + 6, 1)
        self.vtimes = [v_times, v_days]
        # dates
        self.dates = [int(conf_df.yr[0]), int(conf_df.mth[0]),
                      int(conf_df.md[0]), int(conf_df.TT[0])]
        self.track = [float(conf_df.domain_rad[0]), conf_df.track_data_root[0],
                      conf_df.track_data_metadata[0]]
        # Determine the dimensions of the stamp plot according to no members
        nrows, ncols = find_subplot_dims(len(self.ens_members))
        self.rows_cols = [nrows, ncols]
        print('Loaded configuration settings')
        # add stash codes
        self.stashfile = stashfile

    def stash_vars(self):
        """stash_vars
        Description
        Args:
        Returns
        """
        stash_df = pd.read_csv(self.stashfile + '.csv')
        u_stash = stash_df.u_stash[0]
        v_stash = stash_df.v_stash[0]
        u_constraint = iris.AttributeConstraint(STASH=u_stash)
        v_constraint = iris.AttributeConstraint(STASH=v_stash)
        p_constraint = iris.Constraint(pressure=self.plev)
        return u_constraint, v_constraint, p_constraint

    def windloop(self):
        """windloop
        Description:
            Loops over days and hours from config file plotting WindSpeed stamp
            plots
        Args:
            none
        Returns:
            WindSpeed stamp plots
        """
        v_times, v_days = self.vtimes
        for day in v_days:
            for hour in v_times:
                self.ws_dayhour(day, hour)

    def ws_dayhour(self, day, hour):
        """ws_dayhour
        Description:
            Calls the plotter for WindSpeed for set time and adds approriate
            lables
        Args:
            day (int): day
            hour (int): hour
        Returns:
            WindSpeed stamp plot for certain time
        """
        outtest = self.time_const(day, hour)
        year, mth, mmdd, time = self.dates
        if outtest is None or outtest is False:
            return
        # Create figure
        outfile, tcon, bcon, ltime = outtest
        nrows, ncols = self.rows_cols
        fig, axs = plt.subplots(nrows, ncols, dpi=100,
                                subplot_kw={'projection': ccrs.PlateCarree()})
        for i, ens in enumerate(self.ens_members):
            # Determine figure coordinate
            col_i = i % ncols
            row_i = i / ncols
            if len(axs.shape) > 1:
                ax_ws = axs[row_i][col_i]
            else:
                ax_ws = axs[col_i]
            llab, blab = label_fixer(i, ncols, nrows)
            constraints = [tcon, bcon]
            labels = [llab, blab]
            data_f = (self.data_loc.format(mmdd, time)
                      + self.data_name.format(mmdd, time))
            wspeed_plt = self.plot_ens(ax_ws, [data_f, ens],
                                       constraints, labels)
        # Reduce white space and then make whole figure bigger,
        # keeping the aspect ratio constant.
        plt.gcf().subplots_adjust(hspace=0.025, wspace=0.025, bottom=0.05,
                                  top=0.95, left=0.075, right=0.925)
        width, height = plt.gcf().get_size_inches()
        plt.gcf().set_size_inches(width * 3, height * 3)
        # Format colourbar
        cbar = plt.colorbar(wspeed_plt, ax=fig.get_axes(),
                            orientation='horizontal', extend='both',
                            fraction=0.046, pad=0.09)
        cbar.ax.tick_params(labelsize=18)
        cbar.set_label(str(self.plev) + 'hPa WindSpeed (ms$^{-1}$)', size=18)
        current = fig.gca()
        string1 = ('Initial time: {0}/{1}/{2}, {3:02d}' +
                   'Z').format(str(mmdd)[-2:], str(mmdd)[:2], year, time)
        xy1 = [0.95, 0.95]
        string2 = ('Valid time: {0:02d}/{1:02d}/{2}, {3:02d}' +
                   'Z').format(day, mth, year, hour)
        xy2 = [0.95, 0.925]
        string3 = 'T+{0}Z'.format(ltime)
        xy3 = [0.95, 0.9]
        tct.annotate(current, string1, xy1)
        tct.annotate(current, string2, xy2)
        tct.annotate(current, string3, xy3)
        plt.savefig(outfile)
        plt.close()

    def plot_ens(self, ax_ws, fpat_ens, contraints, labels):
        """plot_ens
        Description:
            Plots each subplot for ensemble wind speed
        Args:
            ax_ws (matplotlib figure axis): the axes to plot to
            fpat (str): file path to data
            ens (int): ensemble memeber
            contraints (tuple): tcon (iris time contraint)
                                bcon (iris location constraint)
            labels (tuple): llab (str or logical), left lable or false
                            blab (str or logical), bottome label or false
        Returns:
            approriate subplot
        """
        d_file = fpat_ens[0] + '{0:02d}.pp'.format(fpat_ens[1])
        # Load the data for this ensemble member at this time
        var_constraints = self.stash_vars()
        uvel = vel_extract(d_file, var_constraints[0], contraints[0],
                           var_constraints[2], contraints[1])
        vvel = vel_extract(d_file, var_constraints[1], contraints[0],
                           var_constraints[2], contraints[1])
        uvel, vvel, w_speed = winds(uvel, vvel, fpat_ens[1])
        wspeed_plt = plot_wspeed(ax_ws, w_speed)
        plot_winds(ax_ws, uvel, vvel, self.model)
        # Add grid lines and ticks etc.
        map_formatter(ax_ws, bottom_label=labels[1], left_label=labels[0],
                      labelsize=20, tick_base_x=2, tick_base_y=2)

        # Annotate to add the ensemble member onto the plot
        annotate('Em{0:02d}'.format(fpat_ens[1]), xy=(0.97, 0.03),
                 xycoords='axes fraction',
                 horizontalalignment='right',
                 verticalalignment='bottom', color='k',
                 backgroundcolor='white', fontsize=20)
        return wspeed_plt

    def time_const(self, day, hour):
        """time_const
        Description:
            Checks and generates cube constraints
        Args:
            year (int): year
            month (int): month
            day (int): day
            hour (int): hour
        Returns:
            outfile (str): name of file to be produced
            tcon (iris constraint): time constraint
            bcon (iris constraint): location constraint
            ltime (int): lead time (if false returns nothing)
        """
        year, mth, mmdd, time = self.dates
        ltime = checker(day, hour, self.times)
        if ltime is False:
            print('Skipping: dd = {0}, hour = {1}'.format(day, hour))
            return ltime
        tcon = iris.Constraint(time=iris.time.PartialDateTime(year=year,
                                                              month=mth,
                                                              day=day,
                                                              hour=hour))
        outfile = self.outfile_loc + 'wind/' + self.outfile_name[0]
        outfile = outfile.format(mmdd, time)
        outfile = (outfile + self.model + '_' + str(self.plev)
                   + 'hPa_{0:02d}_{1:02d}Z.png').format(day, hour)
        # Fix domain according to position of ensemble member 0,
        # this assumes the position is similar in ensemble members
        mmll = find_centre_3h([mmdd, day, hour, 0, self.model], self.track)
        if mmll is False:
            print('Skipping: dd = {0}, hour = {1}'.format(day, hour))
            return mmll
        bcon = tct.box_constraint(mmll[0], mmll[1], mmll[2], mmll[3])
        return outfile, tcon, bcon, ltime
