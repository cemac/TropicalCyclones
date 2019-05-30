# -*- coding: utf-8 -*-
"""olr_stamps

.. module:: Outgoing Longwave Radiation Stamp Plotter
    :platform: Unix
    :synopis: Wind Diagnostics and OLR

.. moduleauthor: John Ashcroft, Will Torgerson, Helen Burns CEMAC (UoL)

.. date: May 2019.

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
import cPickle as pickle
import iris
import iris.analysis.calculus
import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as image
import matplotlib.ticker as mticker
import toolkit as tct
mpl.pyplot.switch_backend('Agg')


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
        self.dates = ['03', '09', '2017', 14]
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
        fig = self.loop_olr(self.dates[-1])
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
        plt.rcParams['xtick.labelsize'] = 14
        plt.rcParams['ytick.labelsize'] = 14
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
        time = self.dates[-1]
        olr_ax.set_extent([latlons[1][time] - 5, latlons[1][time] + 5,
                           latlons[0][time] - 5, latlons[0][time] + 5])
        vvcontour = olr_ax.contourf(x, y, dataarray, levels=self.levelsvv,
                                    cmap=plt.cm.get_cmap('binary'),
                                    extend='both')
        cbar_ax = fig.add_axes([0.92, 0.15, 0.05, 0.7])
        cbar = fig.colorbar(vvcontour, cax=cbar_ax)
        cbar.ax.set_title(r'$\mathregular{W/m^{2}}$')
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
        fig.suptitle('Outgoing long wave radiation', fontsize=18)
        day, mnth, year, time = self.dates
        string1 = ('Valid: {0}/{1}/{2} (T+{3}hr)').format(day, mnth, year,
                                                          str(time))
        xy1 = [0.95, 0.95]
        tct.annotate(fig.gca(), string1, xy1)
        plt.savefig(('plots/olr/olr_{0}{1}{2}_{3}Z.png').format(day, mnth,
                                                                year,
                                                                str(time)))
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
