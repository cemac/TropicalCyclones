"""UNIT TEST for windspeed.py
    produce 1 plot,
    must select sensible year, month, day , time

    Notes:
        Very slow
        Original test didn't work for all cases...
    Speed1: 24 mins
    Updates:
    Speed2:
"""

import windspeed

c = windspeed.WindSpeed()
c.dayhour(2014, 12, 04, 18, c.init_day, c.final_day, c.init_time, c.final_time)
