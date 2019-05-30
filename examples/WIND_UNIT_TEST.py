"""UNIT TEST for Winds plot
"""

import numpy as np
import sys
sys.path.append('../TropicalCyclones_py2')
import windspeed_stamps as winds

c = winds.WindSpeedPlotter()
day = 4  # try for non skipped
hour = 12
c.ws_dayhour(day, hour)