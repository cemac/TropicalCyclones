# -*- coding: utf-8 -*-
"""Diags

.. module:: windspeed test
    :platform: Unix
    :synopis: Produces a simple stamp plot of windspeed for an ensemble of size
              N

.. moduleauthor: John Ashcroft, CEMAC (UoL) February 2019.

.. description: This module was developed by CEMAC as part of the WCSSP
                Project. Intial script improvements

   :copyright: © 2019 University of Leeds.
   :license: BSD-2 Clause.

Example:
    To use::

Memebers:

.. CEMAC_stomtracking:
   https://github.com/cemac/TropicalCyclones
"""

import numpy as np
import iris
import ast
import matplotlib as mpl
import cartopy.crs as ccrs
import pandas as pd
import cPickle as pickle
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import sys
import iris.analysis.calculus
import toolkit as tct
# University System python may be broken
# If some one insists on using it...
backend = mpl.get_backend()
if backend == 'Qt4Agg' and sys.version_info[0] == 2:
    # Fix the backend
    print('swapping to Agg Backend')
    print('Please consider using anaconda')
    mpl.use('Agg')
# DO NOT MOVE ABOVE BACKEND FIX
import matplotlib.pyplot as plt # KEEP ME HERE!!!
import matplotlib.image as image
import matplotlib.ticker as mticker


class DiagPlotter(object):
    '''Diagnostics Plotter
        Initialised extra varible and provides standard plots and
        plotting tools

    Members:

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
        Storm = conf_df.Storm[0]
        suitno = conf_df.suitno[0]
        run = conf_df.run[0]
        mmdd_hr = conf_df.mmdd_hr[0]
        var = conf_df.vars[0]
        self.data_loc = (self.data_root + '/' + Storm + '/' + suitno + '/' +
                         run + '/' + mmdd_hr + '/' + var)
        self.data_name = conf_df.data_name[0]
        self.outfile_loc = conf_df.outfile_loc[0]
        self.outfile_name = conf_df.outfile_name.values
        self.ens_members = np.arange(int(conf_df.ens_members[0]))
        self.plev = int(conf_df.plev[0])
        self.md = int(conf_df.md[0])
        self.TT = int(conf_df.TT[0])
        self.model = conf_df.model[0]
        self.v_times = ast.literal_eval(conf_df.v_times[0])
        self.init_day = int(conf_df.init_day[0])
        self.init_time = int(conf_df.init_time[0])
        self.final_day = int(conf_df.final_day[0])
        self.final_time = int(conf_df.final_time[0])
        self.v_days = np.arange(self.init_day, self.init_day + 6, 1)
        self.yr = int(conf_df.yr[0])
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
        self.df = self.data_loc.format(self.md, self.TT) + self.data_name.format(self.md, self.TT)
        print('Loaded stash codes')

        # Plotting configuration
        self.plot_df = pd.read_csv(plotfile + '.csv')
        print('Loaded plot settings file')
        # Read configuration file to import settings
        self.ofile = 'plots/hovs/{0:04d}_{1:02d}Zhovmoller_4p4_{0:04d}_{1:02d}Z_em{2:02d}.png'
        # imfile, olrhome, t select root should be added to config
        imfile = '20170903IRMA.png'
        self.olrhome = '/nfs/a319/ee16wst/OLR/'
        self.timeselect = 14
        self.imfile = self.olrhome + imfile
        self.root = "/nfs/a299/TCs/maria/MARIA_09{1:02}_{2:02}Z_em{0:02}_pb.pp"
        self.levelsvv = ((np.arange(25)) * 10) + 100

    def windoop(self):
        """loop
        Args:
        Returns:
        """
        for dd in self.v_days:
            for hr in self.v_times:
                self.dayhour(self.yr, self.mth, dd, hr, self.init_day,
                             self.final_day, self.init_time, self.final_time)

    def olrloop(self):
        """loop
        Args:
        Returns:
        """
        fig = self.plots_loop(self.timeselect)
        self.finplot(fig)

    def hovloop(self):
        """loop
        Args:
        Returns:
        """
        for md in [self.md]:
            for TT in [self.TT]:
                for em in self.ens_members:
                    self.hovplotter(md, TT, em)

    def dayhour(self, yr, mm, dd, hr, d0, dN, t0, tN):
        """dayhour
        Args:
        Returns:
        """

        outtest = self.t_const(yr, mm, dd, hr, d0, dN, t0, tN)
        if outtest is None:
            return
        # Create figure
        outfile, tcon, bcon, ltime = outtest
        fig, axs = plt.subplots(self.nrows, self.ncols, dpi=100, subplot_kw={
                                'projection': ccrs.PlateCarree()})
        for i, em in enumerate(self.ens_members):
            # Determine figure coordinate
            ic = i % self.ncols
            ir = i / self.ncols
            if len(axs.shape) > 1:
                ax = axs[ir][ic]
            else:
                ax = axs[ic]
            llab, blab = tct.label_fixer(i, self.ncols, self.nrows)
            wspeed_plt = self.plot_i(ax, self.df, em, tcon, bcon, llab, blab)
        # Reduce white space and then make whole figure bigger,
        # keeping the aspect ratio constant.
        plt.gcf().subplots_adjust(hspace=0.025, wspace=0.025, bottom=0.05,
                                  top=0.95, left=0.075, right=0.925)
        w, h = plt.gcf().get_size_inches()
        plt.gcf().set_size_inches(w * 3, h * 3)
        ax = fig.get_axes()
        # Format colourbar
        cbar = plt.colorbar(wspeed_plt, ax=fig.get_axes(),
                            orientation='horizontal', extend='both',
                            fraction=0.046, pad=0.09)
        cbar.ax.tick_params(labelsize=18)
        cbar.set_label(str(self.plev) + 'hPa WindSpeed (ms$^{-1}$)', size=18)
        axs = fig.gca()
        string1 = ('Initial time: {0}/{1}/{2}, {3:02d}'
                   + 'Z').format(str(self.md)[-2:], str(self.md)[:2],
                                 self.yr, self.TT)
        xy1 = [0.95, 0.95]
        string2 = ('Valid time: {0:02d}/{1:02d}/{2}, {3:02d}'
                   + 'Z').format(dd, self.mth, self.yr, hr)
        xy2 = [0.95, 0.925]
        string3 = 'T+{0}Z'.format(ltime)
        xy3 = [0.95, 0.9]
        tct.annotate(axs, string1, xy1)
        tct.annotate(axs, string2, xy2)
        tct.annotate(axs, string3, xy3)
        plt.savefig(outfile)
        plt.close()

    def plot_i(self, ax, fpat, em, tcon, bcon, llab, blab):
        """plot_i
        Args:
        Returns:
        """
        df = fpat + '{0:02d}.pp'.format(em)
        # Load the data for this ensemble member at this time
        u = tct.uv(df, self.u_constraint, tcon, self.p_constraint, bcon)
        v = tct.uv(df, self.v_constraint, tcon, self.p_constraint, bcon)
        u, v, ws = tct.winds(u, v, em)
        wspeed_plt = tct.plot_wspeed(ax, ws)
        tct.plot_winds(ax, u, v, self.model)
        # Add grid lines and ticks etc.
        tct.map_formatter(ax, bottom_label=blab, left_label=llab, labelsize=20,
                          tick_base_x=2, tick_base_y=2)

        # Annotate to add the ensemble member onto the plot
        ax.annotate('{0:02d}'.format(em), xy=(0.97, 0.03),
                    xycoords='axes fraction',
                    horizontalalignment='right',
                    verticalalignment='bottom', color='k',
                    backgroundcolor='white', fontsize=20)
        return wspeed_plt

    def t_const(self, yr, mm, dd, hr, d0, dN, t0, tN):
        """t_const
        Args:
        Returns:
        """
        ltime = tct.checker(dd, hr, d0, dN, t0, tN)
        if ltime is False:
            print('Skipping: dd = {0}, hr = {1}'.format(dd, hr))
            return
        tcon = iris.Constraint(time=iris.time.PartialDateTime(year=yr,
                                                              month=mm, day=dd,
                                                              hour=hr))
        out = self.outfile_loc + 'wind/' + self.outfile_name[0] + '.png'
        outfile = out.format(dd, hr)
        # Fix domain according to position of ensemble member 0,
        # this assumes the position is similar in ensemble members
        mmll = tct.find_centre_3h(self.md, self.TT, 0, dd, hr, self.model,
                                  self.froot, self.fpat, self.domain_rad)
        if mmll is False:
            print('Skipping: dd = {0}, hr = {1}'.format(dd, hr))
            return
        bcon = tct.box_constraint(mmll[0], mmll[1], mmll[2], mmll[3])
        return outfile, tcon, bcon, ltime

    def hovplotter(self, md, TT, em):
        ofile = self.ofile
        fpath = self.data_loc.format(md, TT) + self.data_name.format(md, TT)
        outfile = ofile.format(md, TT, em)
        filepath = self.froot + \
                        '_4p4_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(
                            md, TT, em)
        track_data = np.load(filepath)
        lats = track_data['lats']
        lons = track_data['lons']
        [y0, x0] = [lats, lons]  # centre of storm
        vtan, vrad = tct.load_ens_members(em, fpath, x0, y0, self.u_constraint,
                                          self.v_constraint, self.p_constraint)
        tct.plot_hovmoller(vtan, vrad, outfile, em)

    def plots_loop(self, time):
        """loop
        Args:
        Returns:
        """
        fig = plt.figure(figsize=(15, 12))
        for i in range(18):
            filename = self.olrhome + "olralle{0}.pp".format(i)
            ir = iris.load(filename)[0][0]
            with open(self.olrhome + 'late{0:02}.pkl'.format(i), "rb") as fp:
                latt = np.asarray(pickle.load(fp))
            with open(self.olrhome + 'lone{0:02}.pkl'.format(i), "rb") as fp:
                lonn = np.asarray(pickle.load(fp))
            minlon = float(lonn[time] - 5)
            maxlon = float(lonn[time] + 5)
            minlat = float(latt[time] - 5)
            maxlat = float(latt[time] + 5)
            if i == 0:
                ir_temp = tct.extracter(ir, minlon, maxlon, minlat, maxlat)
            ir = tct.extracter(ir, minlon, maxlon, minlat, maxlat)
            n = i + 1
            newir = self.plotORL(ir, fig, n, latt, lonn, ir_temp)
        return fig

    def plotORL(self, ir, fig, n, latt, lonn, irtemp):
        dataarray = np.zeros((750, 1000))
        dataarray = ir.data
        x = ir.coord('longitude').points
        y = ir.coord('latitude').points
        X, Y = np.meshgrid(x, y)
        ax = fig.add_subplot(4, 5, n + 1, projection=ccrs.PlateCarree())
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=1, color='k', linestyle=':')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlocator = mticker.MultipleLocator(base=2)
        gl.ylocator = mticker.MultipleLocator(base=2)
        gl.xlabel_style = {'size': 8}
        gl.ylabel_style = {'size': 8}
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        ax.coastlines(resolution='10m', color='k', linewidth=1)
        ax.annotate('Em {0}'.format(n - 1), xy=(0.97, 0.03),
                    xycoords='axes fraction', horizontalalignment='right',
                    verticalalignment='bottom', color='k',
                    backgroundcolor='white', fontsize=12)
        ts = self.timeselect
        ax.set_extent([lonn[ts] - 5, lonn[ts] + 5, latt[ts] - 5, latt[ts] + 5])
        vvcontour = ax.contourf(X, Y, dataarray, levels=self.levelsvv,
                                cmap=plt.cm.get_cmap('binary'), extend='both')
        cbar_ax = fig.add_axes([0.92, 0.15, 0.05, 0.7])
        cbar = fig.colorbar(vvcontour, cax=cbar_ax)
        cbar.ax.set_title(r'$\mathregular{W/m^{2}}$')
        newir = ir.regrid(irtemp, iris.analysis.Linear())
        return newir

    def finplot(self, fig):
        ax = fig.add_subplot(4, 5, 1)
        ax.imshow(image.imread(self.imfile))
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.annotate('IR', xy=(0.97, 0.03), xycoords='axes fraction',
                    horizontalalignment='right', verticalalignment='bottom',
                    color='k', backgroundcolor='white', fontsize=12)
        plt.text(x=0.5, y=0.96, s="Outgoing long wave radiation",
                 fontsize=18, ha="center", transform=fig.transFigure)
        plt.text(x=0.5, y=0.912, s="Valid: 03/09/17 14Z (T+14h)",
                 fontsize=12, ha="center", transform=fig.transFigure)
        plt.savefig('test.png')
        plt.close()
