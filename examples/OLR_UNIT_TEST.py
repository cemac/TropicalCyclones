"""UNIT TEST for OLR plot
"""
import numpy as np
import sys
sys.path.append('../TropicalCyclones_py2')
import olr_stams as olr

o = olr.OLR()
olr.olrloop()
