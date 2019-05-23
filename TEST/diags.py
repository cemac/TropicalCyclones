# -*- coding: utf-8 -*-
"""Diags

.. module:: Diagnostics Plotter
    :platform: Unix
    :synopis: Wind Diagnostics and OLR

.. moduleauthor: John Ashcroft, Sam Hardy, Will CEMAC (UoL) February 2019.

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


def set_ofile():
    """set_ofile
    """
    ofile = ('plots/hovs/{0:04d}_{1:02d}Zhovmoller_4p4_{0:04d}_{1:02d}Z_' +
             'em{2:02d}.png')
    return ofile


class OLR(object):
    '''OLR
        Initialised extra varible and provides standard plots and
        plotting tools

    Members:
        olrloop: loop over OLR ensemble members
        loop_olr:
        plot_olr:
        fin_olr_plot:
    '''
    def __init__(self):
        """OLR
        """
        imfile = '20170903IRMA.png'
        self.olrhome = '/nfs/a319/ee16wst/OLR/'
        self.imfile = self.olrhome + imfile
        self.timeselect = 14
        self.levelsvv = ((np.arange(25)) * 10) + 100

    def olrloop(self):
        """olrloop
        Description:
            Execute loop over ensemble plotting of Outgoing longwave radiation
            and add finishing touches
        Args:
            none
        Returns:
            OLR plot for set time accross ensembles.
        """
        fig = self.loop_olr(self.timeselect)
        self.fin_olr_plot(fig)

    def plot_olr(self, cubes, fig, i, latlons):
        """plot_olr
        Description:
            Outgoing Longwave Radiation Plotter.
        Args:
            cubes (tuple of iris cubes):
                                        cube_i constrained iris cube (emN)
                                        cube_temp constrained iris cube (em0)
            fig (matplotlib figure):
            n (int): iterable (plot/ ensemble number)
            latlons (tuple): latt (list) latitude list
                             lonn (list): longitude list
            olr (object): OLR class containing olr info

        Returns:
            fig (matplotlib figure): plot of OLR for set time for all members
        """
        dataarray = np.zeros((750, 1000))
        dataarray = cubes[0].data
        x = cubes[0].coord('longitude').points
        y = cubes[0].coord('latitude').points
        x, y = np.meshgrid(x, y)
        olr_ax = fig.add_subplot(4, 5, i + 1, projection=ccrs.PlateCarree())
        grl = olr_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                               linewidth=1, color='k', linestyle=':')
        grl.xlabels_top = False
        grl.ylabels_right = False
        grl.xlocator = mticker.MultipleLocator(base=2)
        grl.ylocator = mticker.MultipleLocator(base=2)
        grl.xlabel_style = {'size': 8}
        grl.ylabel_style = {'size': 8}
        grl.xformatter = LONGITUDE_FORMATTER
        grl.yformatter = LATITUDE_FORMATTER
        olr_ax.coastlines(resolution='10m', color='k', linewidth=1)
        olr_ax.annotate('Em {0}'.format(i - 1), xy=(0.97, 0.03),
                        xycoords='axes fraction', horizontalalignment='right',
                        verticalalignment='bottom', color='k',
                        backgroundcolor='white', fontsize=12)
        time = self.timeselect
        olr_ax.set_extent([latlons[1][time] - 5, latlons[1][time] + 5,
                           latlons[0][time] - 5, latlons[0][time] + 5])
        vvcontour = olr_ax.contourf(x, y, dataarray, levels=self.levelsvv,
                                    cmap=plt.cm.get_cmap('binary'),
                                    extend='both')
        cbar_ax = fig.add_axes([0.92, 0.15, 0.05, 0.7])
        cbar = fig.colorbar(vvcontour, cax=cbar_ax)
        cbar.olr_ax.set_title(r'$\mathregular{W/m^{2}}$')
        new_cube = cubes[0].regrid(cubes[1], iris.analysis.Linear())
        return new_cube

    def fin_olr_plot(self, fig):
        """fin_olr_plot
        Description:
            Add observation to OLR plot
        Args:
            fig (matplotlib figure): figure axes to add to
        Returns:
            matplotlib figure OLR and saved
        """
        ax_sat = fig.add_subplot(4, 5, 1)
        ax_sat.imshow(image.imread(self.imfile))
        ax_sat.get_xaxis().set_visible(False)
        ax_sat.get_yaxis().set_visible(False)
        ax_sat.annotate('IR', xy=(0.97, 0.03), xycoords='axes fraction',
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        color='k', backgroundcolor='white', fontsize=12)
        plt.text(x=0.5, y=0.96, s="Outgoing long wave radiation",
                 fontsize=18, ha="center", transform=fig.transFigure)
        plt.text(x=0.5, y=0.912, s="Valid: 03/09/17 14Z (T+14h)",
                 fontsize=12, ha="center", transform=fig.transFigure)
        plt.savefig('OLR.png')
        plt.close()

    def loop_olr(self, time):
        """loop_olr
        Description:
            Loops over ensemble members to produce a stamp plot of outgoing
            longwave radiation. Calls data extracter and plotter.
        Args:
            time (int): selected time
        Returns:
            fig (matplotlib figure): plot of OLR for set time for all members
        """
        fig = plt.figure(figsize=(15, 12))
        for i in range(18):
            filename = self.olrhome + "olralle{0}.pp".format(i)
            cube_i = iris.load(filename)[0][0]
            with open(self.olrhome + 'late{0:02}.pkl'.format(i), "rb") as pkl:
                latt = np.asarray(pickle.load(pkl))
            with open(self.olrhome + 'lone{0:02}.pkl'.format(i), "rb") as pkl:
                lonn = np.asarray(pickle.load(pkl))
            minlon = float(lonn[time] - 5)
            maxlon = float(lonn[time] + 5)
            minlat = float(latt[time] - 5)
            maxlat = float(latt[time] + 5)
            if i == 0:
                cube_temp = tct.extracter(cube_i, minlon, maxlon, minlat,
                                          maxlat)
            cube_i = tct.extracter(cube_i, minlon, maxlon, minlat, maxlat)
            self.plot_olr([cube_i, cube_temp], fig, i + 1, [latt, lonn])
        return fig


class DiagPlotter(object):
    '''Diagnostics Plotter
        Initialised extra varible and provides standard plots and
        plotting tools

    Members:
        windloop: loop over windspeed plots (ws_dayhour)
        hovloop: loop over ensemble members at give time
        ws_dayhour: loop over the day and hour and call stamp plotter
        plot_ens: plot WindSpeed plot for each ensemble member
        time_const: contrain for time and location
        hovplotter: call hovoller plotter for time and memeber
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
        self.data_loc = (conf_df.data_root[0] + '/' + conf_df.Storm[0] + '/' +
                         conf_df.suitno[0] + '/' + conf_df.run[0] + '/' +
                         conf_df.mmdd_hr[0] + '/' + conf_df.vars[0])
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
        v_times = ast.literal_eval(conf_df.v_times[0])
        v_days = np.arange(init_day, init_day + 6, 1)
        self.vtimes = [v_times, v_days]
        # dates
        self.dates = [int(conf_df.year[0]), int(conf_df.mth[0]),
                      int(conf_df.md[0]), int(conf_df.TT[0])]
        self.track = [float(conf_df.domain_rad[0]), conf_df.track_data_root[0],
                      conf_df.track_data_metadata[0]]
        # Determine the dimensions of the stamp plot according to no members
        nrows, ncols = tct.find_subplot_dims(len(self.ens_members))
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

    def hovloop(self):
        """hovloop
        Description:
            Single plots for each ensemble member of tangential, radial and
            vertical wind speed.
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
        if outtest is None:
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
            llab, blab = tct.label_fixer(i, ncols, nrows)
            constraints = [tcon, bcon]
            labels = [llab, blab]
            data_f = (self.data_loc.format(mmdd, time) +
                      self.data_name.format(mmdd, time))
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
        string1 = ('Initial time: {0}/{1}/{2}, {3:02d}'
                   + 'Z').format(str(mmdd)[-2:], str(mmdd)[:2], year, time)
        xy1 = [0.95, 0.95]
        string2 = ('Valid time: {0:02d}/{1:02d}/{2}, {3:02d}'
                   + 'Z').format(day, mth, year, hour)
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
        uvel = tct.vel_extract(d_file, var_constraints[0], contraints[0],
                               var_constraints[2], contraints[1])
        vvel = tct.vel_extract(d_file, var_constraints[1], contraints[0],
                               var_constraints[2], contraints[1])
        uvel, vvel, w_speed = tct.winds(uvel, vvel, fpat_ens[1])
        wspeed_plt = tct.plot_wspeed(ax_ws, w_speed)
        tct.plot_winds(ax_ws, uvel, vvel, self.model)
        # Add grid lines and ticks etc.
        tct.map_formatter(ax_ws, bottom_label=labels[1], left_label=labels[0],
                          labelsize=20, tick_base_x=2, tick_base_y=2)

        # Annotate to add the ensemble member onto the plot
        ax_ws.annotate('{0:02d}'.format(fpat_ens[1]), xy=(0.97, 0.03),
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
        ltime = tct.checker(day, hour, self.times)
        if ltime is False:
            print('Skipping: dd = {0}, hour = {1}'.format(day, hour))
            return ltime
        tcon = iris.Constraint(time=iris.time.PartialDateTime(year=year,
                                                              month=mth,
                                                              day=day,
                                                              hour=hour))
        out = self.outfile_loc + 'wind/' + self.outfile_name[0] + '.png'
        outfile = out.format(day, hour)
        # Fix domain according to position of ensemble member 0,
        # this assumes the position is similar in ensemble members
        mmll = tct.find_centre_3h([mmdd, time, 0, day, hour, self.model],
                                  self.v_times, self.track)
        if mmll is False:
            print('Skipping: dd = {0}, hour = {1}'.format(day, hour))
            return mmll
        bcon = tct.box_constraint(mmll[0], mmll[1], mmll[2], mmll[3])
        return outfile, tcon, bcon, ltime

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
        filepath = self.track[0] + \
            '_4p4_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(
                mmdd, hour, ens)
        track_data = np.load(filepath)
        lats = track_data['lats']
        lons = track_data['lons']
        [y_0, x_0] = [lats, lons]  # centre of storm
        contraints = self.stash_vars()
        vtan, vrad = tct.load_ens_members(ens, fpath, x_0, y_0, constraints)
        tct.plot_hovmoller(vtan, vrad, outfile, ens)
