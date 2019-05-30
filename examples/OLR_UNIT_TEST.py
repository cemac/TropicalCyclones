"""UNIT TEST for OLR plot
"""
import numpy as np
import sys
sys.path.append('../TropicalCyclones_py2')
import olr_stamps as olr

o = olr.OLR()
o.olrloop()
