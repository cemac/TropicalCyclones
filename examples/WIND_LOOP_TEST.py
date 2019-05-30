"""UNIT TEST for Winds plot
"""

import numpy as np
import sys
sys.path.append('../TropicalCyclones_py2')
import windspeed_stamps as winds

c = winds.WindSpeedPlotter()
c.windloop()
