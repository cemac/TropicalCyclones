# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt  # KEEP ME HERE!!!
"""hovmoller plots

.. module:: hov test
    :platform: Unix
    :synopis: Produces a simple hovmoller plot

.. moduleauthor: John Ashcroft, Sam Hardy, Helen Burns CEMAC (UoL) May 2019.

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
import ast
import matplotlib as mpl
import matplotlib.pyplot as plt
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


class hovmoller(object):
    '''hovmoller


    Members:

    '''

    def __init__(self, configfile='configfile'):
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
        self.mmdd_hr = conf_df.mmdd_hr[0]
        var = conf_df.vars[0]
        self.data_loc = (self.data_root + '/' + Storm + '/' + suitno + '/' +
                         run + '/' + self.mmdd_hr + '/' + var)
        self.data_name = conf_df.data_name[0]
        self.outfile_loc = conf_df.outfile_loc[0]
        self.outfile_name = conf_df.outfile_name .values
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
        self.froot = conf_df.track_data_root[0]
        self.fpat = conf_df.track_data_metadata[0]
        print('Loaded configuration settings')

    def hovloop(self):
        """loop
        Args:
        Returns:
        """
        ofile = 'plots/hovs/{0:04d}_{1:02d}Zhovmoller_4p4_{0:04d}_{1:02d}Z_em{2:02d}.png'
        for md in [self.md]:
            for TT in [self.TT]:
                fpath = self.data_loc.format(
                    md, TT) + self.data_name.format(md, TT)
                for em in self.ens_members:
                    outfile = ofile.format(md, TT, em)
                    filepath = self.froot + \
                        '_4p4_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(
                            md, TT, em)
                    track_data = np.load(filepath)
                    lats = track_data['lats']
                    lons = track_data['lons']
                    [y0, x0] = [lats, lons]  # centre of storm
                    df = fpath + '{0:02d}.pp'.format(em)
                    v_azi_list = tct.load_ens_members(em, fpath, x0, y0)
                    tct.plot_hovmoller(v_azi_list, outfile)
