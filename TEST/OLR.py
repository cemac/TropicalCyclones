# -*- coding: utf-8 -*-
"""ORL

.. module:: ORL
    :platform: Unix
    :synopis: Produces a simple stamp plot of ORL for an ensemble of size
              N

.. moduleauthor: Will Torgerson, Helen Burns CEMAC (UoL) February 2019.

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
import cPickle as pickle
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import pandas as pd
import ast
import iris.analysis.calculus
import toolkit as tct
import cartopy.crs as ccrs
import numpy as np
import iris
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
# MUSE IMPORT AFTER - ignore linter
import matplotlib.image as image
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt


class OLR(object):
    '''OLR

    Members:

    '''

    def __init__(self, timeselect, imfile='20170903IRMA.png'):
        '''
        Args:
            timeselect (int): number of hours after run
            stashfile (string): filepath sat image
        '''
        # Read configuration file to import settings
        self.olrhome = '/nfs/a319/ee16wst/OLR/'
        self.timeselect = 14
        self.imfile = self.olrhome + imfile
        self.root = "/nfs/a299/TCs/maria/MARIA_09{1:02}_{2:02}Z_em{0:02}_pb.pp"
        self.levelsvv = ((np.arange(25)) * 10) + 100

    def loop(self):
        """loop
        Args:
        Returns:
        """
        newir, fig = self.plots_loop(self.timeselect)
        self.finplot(newir, fig)

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
            x = ir.coord('longitude').points
            y = ir.coord('latitude').points
            X, Y = np.meshgrid(x, y)
            n = i + 1
            newir = self.plotORL(ir, fig, n, latt, lonn, ir_temp)
        return newir, fig

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

    def finplot(self, newir, fig):
        x2 = newir.coord('longitude').points
        y2 = newir.coord('latitude').points
        X2, Y2 = np.meshgrid(x2, y2)
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
