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
import matplotlib
from matplotlib.colors import from_levels_and_colors
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.cm as mpl_cm
