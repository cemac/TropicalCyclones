"""UNIT TEST for Hovmoller plot
"""

import toolkit as tct
import numpy as np

v_azi = np.load(vazi.npy, v_azi)
vrt = np.load(vrt.npy, vrt)
vrad = np.load(vrad.npy, vrad)
outfile = ''
ens = 01
tct.plot_hovmoller(v_azi, vrad, vrt, outfile, ens)
