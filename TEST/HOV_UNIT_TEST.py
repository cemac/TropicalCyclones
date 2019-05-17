"""UNIT TEST for windspeed.py
    produce 1 plot,
    must select sensible year, month, day , time

    Notes:
        Very slow
    Speed1:
    Updates:
    Speed2:
"""

import hov

c = hov.hovmoller()
ofile = 'plots/hovs/{0:04d}_{1:02d}Zhovmoller_4p4_{0:04d}_{1:02d}Z_em{2:02d}.png'
c.hovplotter(1203, 12, 02, ofile)
# 01 - original
# 02 - adding Radial
