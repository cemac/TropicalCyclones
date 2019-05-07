"""UNIT TEST for windspeed.py
    produce 1 plot,
    must select sensible year, month, day , time

    Notes:
        Very slow
    Speed1:
    Updates:
    Speed2:
"""

import OLR

c = OLR.OLR()
c.dayhour(2014, 12, 04, 18, c.init_day, c.final_day, c.init_time, c.final_time)
