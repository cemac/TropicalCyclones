"""UNIT TEST for Hovmoller plot
"""

import numpy as np
import sys
sys.path.append('../TropicalCyclones_py2')
import hovmoller as hov

# Load in save data
v_azi = np.load('data/vazi_EM01.npy')
vrt = np.load('data/vrt_EM01.npy')
vrad = np.load('data/vrad_EM01.npy')
# Set output name and corresponding ensemble number
outfile = 'plots/1203_12Zhovmoller_4p4_1203_12Z_em01.png'
ens = 1
hov.plot_hovmoller(v_azi, vrad, vrt, outfile, ens)
