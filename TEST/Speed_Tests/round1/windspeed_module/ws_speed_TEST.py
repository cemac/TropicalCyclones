# -*- coding: utf-8 -*-
"""Speed checker

.. module:: Python Speed test
    :platform: Unix
    :synopis: Perform all windspeed plots

.. moduleauthor: CEMAC (UoL) February 2019.

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

import time
import windspeed

print('starting')
start_time = time.time()
c = windspeed.WindSpeed()
c.loop()
print("--- %s seconds ---" % (time.time() - start_time))
print('ending')
