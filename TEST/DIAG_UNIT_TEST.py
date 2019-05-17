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
print('IMPORT')
print('time for import')
print("--- %s seconds ---" % (time.time() - start_time))
print('Class initialisation')
start_time = time.time()
c = diags.DiagPlotter()
print('Load Class:')
print('time for Class initialisation')
print("--- %s seconds ---" % (time.time() - start_time))
print('OLR')
start_time = time.time()
c.olrloop()
print('time for single OLR')
print("--- %s seconds ---" % (time.time() - start_time))
print('WindSpeed:')
start_time = time.time()
day = 4  # try for non skipped
hour = 18
c.WSdayhour(c.yr, c.mth, day, hour, c.init_day, c.final_day,
            c.init_time, c.final_time)
print('WindSpeed time:')
print("--- %s seconds ---" % (time.time() - start_time))
print('Hovs')
em = 01
c.hovplotter(c.md, c.mth, em)
print('Hovs time:')
print("--- %s seconds ---" % (time.time() - start_time))
