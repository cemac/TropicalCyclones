"""UNIT TEST for Hovmoller plot
"""

import numpy as np
import sys
sys.path.append('../TropicalCyclones_py2')
import hovmoller as hov

# Load in save data
v_azi = np.load('vazi_EM01.npy')
vrt = np.load('vrt_EM01.npy')
vrad = np.load('vrad_EM01.npy')
# Set output name and corresponding ensemble number
outfile = 'output/1203_12Zhovmoller_4p4_1203_12Z_em01.png'
ens = 1
hovs.plot_hovmoller(v_azi, vrad, vrt, outfile, ens)
