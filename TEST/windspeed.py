# -*- coding: utf-8 -*-
"""windspeed

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


class WindSpeed(object):
    '''WindSpeed
        Initialised extra varible and provides standard plots and
        plotting tools

    Members:

    '''

    def __init__(self, configfile='configfile', stashfile='stashvars',
                 plotfile='plot_conf'):
        '''
        Args:
            configfile (string): filepath to configuration settings
            stashfile (string): filepath to shashfile codes
        '''
        # Read configuration file to import settings
        conf_df = pd.read_csv(configfile + '.csv')
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

    def loop(self):
        """loop
        Args:
        Returns:
        """
        for dd in self.v_days:
            for hr in self.v_times:
                self.dayhour(self.yr, self.mth, dd, hr, self.init_day,
                             self.final_day, self.init_time, self.final_time)

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
        cbar.set_label('ms$^{-1}$', size=18)
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
            print('skipping')
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
            print('skipping')
            return
        bcon = tct.box_constraint(mmll[0], mmll[1], mmll[2], mmll[3])
        return outfile[0], tcon, bcon, ltime
