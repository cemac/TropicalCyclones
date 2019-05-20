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


class DiagPlotter(object):
    '''Diagnostics Plotter
        Initialised extra varible and provides standard plots and
        plotting tools

    Members:
        windloop: loop over windspeed plots (ws_dayhour)
        olrloop: loop over OLR ensemble members
        hovloop: loop over ensemble members at give time
        ws_dayhour: loop over the day and hour and call stamp plotter
        plot_ens: plot WindSpeed plot for each ensemble member
        time_const: contrain for time and location
        hovplotter: call hovoller plotter for time and memeber
        loop_olr:
        plot_olr:
        fin_olr_plot:
    '''

    def __init__(self, configfile='configfile2', stashfile='stashvars',
                 plotfile='plot_conf'):
        '''
        Args:
            configfile (string): filepath to configuration settings
            stashfile (string): filepath to shashfile codes
        '''
        # Read configuration file to import settings
        conf_df = pd.read_csv(configfile + '.csv')
        conf_df = conf_df.set_index('var:').transpose()
        # Define all constraints and locate data.
        self.data_root = conf_df.data_root[0]
        storm = conf_df.Storm[0]
        suitno = conf_df.suitno[0]
        run = conf_df.run[0]
        mmdd_hr = conf_df.mmdd_hr[0]
        var = conf_df.vars[0]
        self.data_loc = (self.data_root + '/' + storm + '/' + suitno + '/' +
                         run + '/' + mmdd_hr + '/' + var)
        self.data_name = conf_df.data_name[0]
        self.outfile_loc = conf_df.outfile_loc[0]
        self.outfile_name = conf_df.outfile_name.values
        self.ens_members = np.arange(int(conf_df.ens_members[0]))
        self.plev = int(conf_df.plev[0])
        self.mmdd = int(conf_df.md[0])
        self.time = int(conf_df.TT[0])
        self.model = conf_df.model[0]
        self.v_times = ast.literal_eval(conf_df.v_times[0])
        self.init_day = int(conf_df.init_day[0])
        self.init_time = int(conf_df.init_time[0])
        self.final_day = int(conf_df.final_day[0])
        self.final_time = int(conf_df.final_time[0])
        self.v_days = np.arange(self.init_day, self.init_day + 6, 1)
        self.year = int(conf_df.year[0])
        self.mth = int(conf_df.mth[0])
        self.domain_rad = float(conf_df.domain_rad[0])
        print('Loaded configuration settings')
        # add stash codes
        stash_df = pd.read_csv(stashfile + '.csv')
        u_stash = stash_df.u_stash[0]
        v_stash = stash_df.v_stash[0]
        self.slp_stash = stash_df.slp_stash[0]
        self.u_constraint = iris.AttributeConstraint(STASH=u_stash)
        self.v_constraint = iris.AttributeConstraint(STASH=v_stash)
        self.froot = conf_df.track_data_root[0]
        self.fpat = conf_df.track_data_metadata[0]
        # Determine the dimensions of the stamp plot according to no members
        self.n_ems = len(self.ens_members)
        self.nrows, self.ncols = tct.find_subplot_dims(self.n_ems)
        self.p_constraint = iris.Constraint(pressure=self.plev)
        self.data_f = self.data_loc.format(
            self.mmdd, self.time) + self.data_name.format(self.mmdd, self.time)
        print('Loaded stash codes')

        # Plotting configuration
        self.plot_df = pd.read_csv(plotfile + '.csv')
        print('Loaded plot settings file')
        # Read configuration file to import settings
        self.ofile = ('plots/hovs' + '/{0:04d}_{1:02d}Zhovmoller_'
                      + '4p4_{0:04d}_{1:02d}Z_em{2:02d}.png')
        # imfile, olrhome, t select root should be added to config
        imfile = '20170903IRMA.png'
        self.olrhome = '/nfs/a319/ee16wst/OLR/'
        self.timeselect = 14
        self.imfile = self.olrhome + imfile
        self.root = "/nfs/a299/TCs/maria/MARIA_09{1:02}_{2:02}Z_em{0:02}_pb.pp"
        self.levelsvv = ((np.arange(25)) * 10) + 100

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
        for day in self.v_days:
            for hour in self.v_times:
                self.ws_dayhour(self.year, self.mth, day, hour,
                                self.init_day, self.final_day, self.init_time,
                                self.final_time)

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
        for mmdd in [self.mmdd]:
            for hour in [self.time]:
                for ens in self.ens_members:
                    self.hovplotter(mmdd, hour, ens)

    def ws_dayhour(self, year, month, day, hour, day_0, day_n, time_0, time_n):
        """ws_dayhour
        Description:
            Calls the plotter for WindSpeed for set time and adds approriate
            lables
        Args:
            year (int): year
            month (int): month
            day (int): day
            hour (int): hour
            day_0 (int): initial day
            day_n (int): end day
            time_0 (int): initial time
            time_n (int): end time
        Returns:
            WindSpeed stamp plot for certain time
        """

        outtest = self.time_const(
            year, month, day, hour, day_0, day_n, time_0, time_n)
        if outtest is None:
            return
        # Create figure
        outfile, tcon, bcon, ltime = outtest
        fig, axs = plt.subplots(self.nrows, self.ncols, dpi=100,
                                subplot_kw={'projection': ccrs.PlateCarree()})
        for i, ens in enumerate(self.ens_members):
            # Determine figure coordinate
            col_i = i % self.ncols
            row_i = i / self.ncols
            if len(axs.shape) > 1:
                ax_ws = axs[row_i][col_i]
            else:
                ax_ws = axs[col_i]
            llab, blab = tct.label_fixer(i, self.ncols, self.nrows)
            wspeed_plt = self.plot_ens(
                ax_ws, self.data_f, ens, tcon, bcon, llab, blab)
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
                   + 'Z').format(str(self.mmdd)[-2:], str(self.mmdd)[:2],
                                 self.year, self.time)
        xy1 = [0.95, 0.95]
        string2 = ('Valid time: {0:02d}/{1:02d}/{2}, {3:02d}'
                   + 'Z').format(day, self.mth, self.year, hour)
        xy2 = [0.95, 0.925]
        string3 = 'T+{0}Z'.format(ltime)
        xy3 = [0.95, 0.9]
        tct.annotate(current, string1, xy1)
        tct.annotate(current, string2, xy2)
        tct.annotate(current, string3, xy3)
        plt.savefig(outfile)
        plt.close()

    def plot_ens(self, ax_ws, fpat, ens, tcon, bcon, llab, blab):
        """plot_ens
        Description:
            Plots each subplot for ensemble wind speed
        Args:
            ax_ws (matplotlib figure axis): the axes to plot to
            fpat (str): file path to data
            ens (int): ensemble memeber
            tcon (iris constraint): time contraints
            bcon (iris constraint): location contraints
            llab (str or logical): either left lable or false
            blab (str or logical): either bottome label or false
        Returns:
            approriate subplot
        """
        d_file = fpat + '{0:02d}.pp'.format(ens)
        # Load the data for this ensemble member at this time
        uvel = tct.uv(d_file, self.u_constraint, tcon, self.p_constraint, bcon)
        vvel = tct.uv(d_file, self.v_constraint, tcon, self.p_constraint, bcon)
        uvel, vvel, w_speed = tct.winds(uvel, vvel, ens)
        wspeed_plt = tct.plot_wspeed(ax_ws, w_speed)
        tct.plot_winds(ax_ws, uvel, vvel, self.model)
        # Add grid lines and ticks etc.
        tct.map_formatter(ax_ws, bottom_label=blab, left_label=llab,
                          labelsize=20, tick_base_x=2, tick_base_y=2)

        # Annotate to add the ensemble member onto the plot
        ax_ws.annotate('{0:02d}'.format(ens), xy=(0.97, 0.03),
                       xycoords='axes fraction',
                       horizontalalignment='right',
                       verticalalignment='bottom', color='k',
                       backgroundcolor='white', fontsize=20)
        return wspeed_plt

    def time_const(self, year, month, day, hour, day_0, day_n, time_0, time_n):
        """time_const
        Description:
            Checks and generates cube constraints
        Args:
            year (int): year
            month (int): month
            day (int): day
            hour (int): hour
            day_0 (int): initial day
            day_n (int): end day
            time_0 (int): initial time
            time_n (int): end time
        Returns:
            outfile (str): name of file to be produced
            tcon (iris constraint): time constraint
            bcon (iris constraint): location constraint
            ltime (int): lead time (if false returns nothing)
        """
        ltime = tct.checker(day, hour, day_0, day_n, time_0, time_n)
        if ltime is False:
            print('Skipping: dd = {0}, hour = {1}'.format(day, hour))
            return ltime
        tcon = iris.Constraint(time=iris.time.PartialDateTime(year=year,
                                                              month=month,
                                                              day=day,
                                                              hour=hour))
        out = self.outfile_loc + 'wind/' + self.outfile_name[0] + '.png'
        outfile = out.format(day, hour)
        # Fix domain according to position of ensemble member 0,
        # this assumes the position is similar in ensemble members
        mmll = tct.find_centre_3h(self.mmdd, self.time, 0, day, hour,
                                  self.model, self.froot, self.fpat,
                                  self.domain_rad)
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
        ofile = self.ofile
        fpath = self.data_loc.format(
            mmdd, hour) + self.data_name.format(mmdd, hour)
        outfile = ofile.format(mmdd, hour, ens)
        filepath = self.froot + \
            '_4p4_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(
                mmdd, hour, ens)
        track_data = np.load(filepath)
        lats = track_data['lats']
        lons = track_data['lons']
        [y_0, x_0] = [lats, lons]  # centre of storm
        vtan, vrad = tct.load_ens_members(ens, fpath, x_0, y_0,
                                          self.u_constraint, self.v_constraint,
                                          self.p_constraint)
        tct.plot_hovmoller(vtan, vrad, outfile, ens)

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
            plot_no = i + 1
            self.plot_olr(cube_i, fig, plot_no, latt, lonn, cube_temp)
        return fig

    def plot_olr(self, cube_i, fig, i, latt, lonn, cube_temp):
        """plot_olr
        Description:
            Outgoing Longwave Radiation Plotter.
        Args:
            cube_i (iris cube):
            fig (matplotlib figure):
            n (int): iterable (plot/ ensemble number)
            latt (list) latitude list
            lonn (list): longitude list
            cube_temp (iris cube): constrained iris cube (em0)
        Returns:
            fig (matplotlib figure): plot of OLR for set time for all members
        """
        dataarray = np.zeros((750, 1000))
        dataarray = cube_i.data
        x_points = cube_i.coord('longitude').points
        y_points = cube_i.coord('latitude').points
        x, y = np.meshgrid(x_points, y_points)
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
        olr_ax.set_extent([lonn[time] - 5, lonn[time] + 5, latt[time] - 5,
                           latt[time] + 5])
        vvcontour = olr_ax.contourf(x, y, dataarray, levels=self.levelsvv,
                                    cmap=plt.cm.get_cmap('binary'),
                                    extend='both')
        cbar_ax = fig.add_axes([0.92, 0.15, 0.05, 0.7])
        cbar = fig.colorbar(vvcontour, cax=cbar_ax)
        cbar.olr_ax.set_title(r'$\mathregular{W/m^{2}}$')
        new_cube = cube_i.regrid(cube_temp, iris.analysis.Linear())
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
