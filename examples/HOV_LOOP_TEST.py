"""UNIT TEST for Hovmoller plot
"""

import numpy as np
import sys
sys.path.append('../TropicalCyclones_py2')
import hovmoller as hov

c = hov.HovPlotter()
ens = 1
c.hovplotter(c.dates[-2], c.dates[-3], ens)
