"""UNIT TEST for Hovmoller plot
"""

import toolkit as tct
import numpy as np

v_azi = np.load('vazi.npy')
vrt = np.load('vrt.npy')
vrad = np.load('vrad.npy')
outfile = 'plots/hovs/1203_12Zhovmoller_4p4_1203_12Z_em01.png'
ens = 1
tct.plot_hovmoller(v_azi, vrad, vrt, outfile, ens)
