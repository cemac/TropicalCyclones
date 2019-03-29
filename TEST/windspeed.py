# -*- coding: utf-8 -*-
"""windspeed

.. module:: S_Plotter
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
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from toolkit import *


class windspeed(object):
    '''StormStats
        Initialised extra varible and provides standard plots and
        plotting tools

    Members:

    '''
    def __init__(self, configfile):

        # Define all constraints and locate data.
        data_loc = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/4p4/1203_12Z/wind_plev/'
        # Should include everything but 'NN.pp' where NN is the ensemble member
        data_name = 'uvplev_4p4_1203_12Z_em'
        outfile_loc = './'
        # time data to be added to the outfile
        outfile_name = 'wspeed_1203_12Z_4p4_850hPa_'

        ens_members = np.arange(12)  # Ensemble members to use (i.e. 12)
        plev = 850  # Must be alevel that is in the data.

        md = 1203
        TT = 12  # Forecast initial times
        model = '4p4'

        v_times = [0, 3, 6, 9, 12, 15, 18, 21]  # Times during day to plot
        init_day = 3
        init_time = 12  # The first plot will be of the time after this time
        final_day = 8
        final_time = 9  # Last time to plot
        v_days = np.arange(init_day, init_day + 6, 1)
        yr = 2014
        mth = 12

        domain_rad = 4.  # Domain is 2*domain_rad x 2*domain_rad

        u_stash = 'm01s15i201'
        v_stash = 'm01s15i202'
        slp_stash = 'm01s16i222'
        u_constraint = iris.AttributeConstraint(STASH=u_stash)
        v_constraint = iris.AttributeConstraint(STASH=v_stash)

        # Determine the dimensions of the stamp plot according to no members
        n_ems = len(ens_members)
        nrows, ncols = find_subplot_dims(n_ems)

        p_constraint = iris.Constraint(pressure=plev)

        for dd in v_days:
            for hr in v_times:
                if dd < init_day or (dd == init_day and hr < init_time + 1):
                    continue
                if dd > final_day or (dd > final_day - 1 and hr > final_time - 1):
                    continue
                lead_time = (dd - init_day) * 24 + (hr - init_time)
                print 'Plotting data for dd = {0}, hr = {1}'.format(dd, hr)
                lead_time = (dd - init_day) * 24 + (hr - init_time)
                time_constraint = iris.Constraint(
                    time=iris.time.PartialDateTime(year=yr, month=mth,
                                                   day=dd, hour=hr))
                outfile = outfile_loc + outfile_name + \
                    '{0:02d}_{1:02d}Z.png'.format(dd, hr)

                # Fix domain according to position of ensemble member 0,
                # this assumes the position is similar in ensemble members
                [cenlat, cenlon] = find_centre_3h(md, TT, 0, dd, hr, model)
                minlat = cenlat - domain_rad
                maxlat = cenlat + domain_rad
                minlon = cenlon - domain_rad
                maxlon = cenlon + domain_rad
                b_constraint = box_constraint(minlat, maxlat, minlon, maxlon)

                # Create figure
                fig, axs = plt.subplots(nrows, ncols, dpi=100, subplot_kw={
                                        'projection': ccrs.PlateCarree()})
                for i, em in enumerate(ens_members):

                    # Determine figure coordinate
                    ic = i % ncols
                    ir = i / ncols
                    if len(axs.shape) > 1:
                        ax = axs[ir, ic]
                    else:
                        ax = axs[ic]

                    # Load the data for this ensemble member at this time
                    df = data_loc + data_name + '{0:02d}.pp'.format(em)
                    with iris.FUTURE.context(cell_datetime_objects=True):
                        u = iris.load_cube(df, u_constraint).extract(
                            time_constraint & p_constraint & b_constraint)
                        v = iris.load_cube(df, v_constraint).extract(
                            time_constraint & p_constraint & b_constraint)
                    ws = (u**2 + v**2)**0.5

                    coord = iris.coords.AuxCoord(em, long_name='ensemble_member')
                    u.add_aux_coord(coord)
                    v.add_aux_coord(coord)
                    ws.add_aux_coord(coord)

                    wspeed_plt = plot_wspeed(ax, ws)
                    plot_winds(ax, u, v, model)

                    # Make sure only the outside axis are labelled
                    if i % ncols != 0:
                        left_label = False
                    else:
                        left_label = True
                    if i < ncols * (nrows - 1):
                        bottom_label = False
                    else:
                        bottom_label = True

                    # Add grid lines and ticks etc.
                    map_formatter(ax, bottom_label=bottom_label,
                                  left_label=left_label, labelsize=20,
                                  tick_base_x=2, tick_base_y=2)

                    # Annotate to add the ensemble member onto the plot
                    ax.annotate('{0:02d}'.format(em), xy=(0.97, 0.03),
                                xycoords='axes fraction',
                                horizontalalignment='right',
                                verticalalignment='bottom', color='k',
                                backgroundcolor='white', fontsize=20)

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

                # Add initial time, valid time etc.
                axs.annotate('Initial time: {0}/{1}/{2}, {3:02d}Z'.format(str(md)[-2:],
                             str(md)[:2], yr, TT), xy=(0.95, 0.95),
                             xycoords='figure fraction',
                             horizontalalignment='right',
                             verticalalignment='top', color='k', fontsize=15)
                axs.annotate('Valid time: {0:02d}/{1:02d}/{2}, {3:02d}Z'.format(dd, mth, yr, hr),
                             xy=(0.95, 0.925), xycoords='figure fraction',
                             horizontalalignment='right', verticalalignment='top',
                             color='k', fontsize=15)
                axs.annotate('T+{0}Z'.format(lead_time), xy=(0.95, 0.9),
                             xycoords='figure fraction',
                             horizontalalignment='right',
                             verticalalignment='top', color='k', fontsize=15)

                plt.savefig(outfile)
                plt.close()
