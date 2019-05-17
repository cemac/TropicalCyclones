"""UNIT TEST for windspeed.py
    produce 1 plot,
    must select sensible year, month, day , time

    Notes:
        Very slow
    Speed1:
    Updates:
    Speed2:
"""
import time
import hov




print('starting')
start_time = time.time()
c = hov.hovmoller()
ofile = 'plots/hovs/{0:04d}_{1:02d}Zhovmoller_4p4_{0:04d}_{1:02d}Z_em{2:02d}.png'
c.hovloop()
print("--- %s seconds ---" % (time.time() - start_time))
print('ending')
