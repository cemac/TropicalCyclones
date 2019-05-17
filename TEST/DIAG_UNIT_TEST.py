"""UNIT TEST for diags.py
    produce 1 plot,
    must select sensible year, month, day , time

    Notes:
        Very slow
    Speed1:
    Updates:
    Speed2:
"""
import time
start_time = time.time()
import diags
print('time for import')
print("--- %s seconds ---" % (time.time() - start_time))
print('Class initialisation')
start_time = time.time()
c = diags.DiagPlotter()
print('time for Class initialisation')
print("--- %s seconds ---" % (time.time() - start_time))
print('OLR')
start_time = time.time()
ofile = 'plots/hovs/{0:04d}_{1:02d}Zhovmoller_4p4_{0:04d}_{1:02d}Z_em{2:02d}.png'
c.olrloop()
print('time for single OLR')
print("--- %s seconds ---" % (time.time() - start_time))
print('WindSpeed')
start_time = time.time()
c.dayhour(2014, 12, 01, 18, c.init_day, c.final_day, c.init_time, c.final_time)
print('WindSpeed time:')
print("--- %s seconds ---" % (time.time() - start_time))
print('Hovs')
c.hovplotter(1203, 12, 02, ofile)
print('Hovs time:')
print("--- %s seconds ---" % (time.time() - start_time))
# 01 - original
# 02 - adding Radial
