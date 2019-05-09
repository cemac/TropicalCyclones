########################################################################
# Selection of plotting tools specific to my needs. Contents below.    #
#                                              John Ashcroft, Oct 2017 #
########################################################################

import cartopy.crs as ccrs
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.cm as mpl_cm

import iris
#import iris.quickplot as qplt
from matplotlib.colors import from_levels_and_colors
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.colors as colors
import cmocean
#import cmocean
#### Contents:   #####
#
## map_formatter()
# Draws gridlines, coastlines, ticks and formats labels on axes. 
# tick_base defines the interval in degrees between points on the grid.
#
## normalise_cmap()
# inputs a max, min, midpoint and spacing (dl). The colourmap is then
# centred on the midpoint. Used in particular for diverging colourmaps
# so that the white colour is equal to the value of 0.
# Returns cmap, norm and levels (i.e. an example of the plot after using
# function: ax.contourf(X,Y,DATA,norm=norm,levels=levels,cmap=cmap)
#
## OOMFormatter()
# Used in scripts to format how scientific numbers are represented. 
# Ensures that everything is neat and tidy with the correct power abover
# the colourbar. For example, see vorticity.py
#
## highResPoints()
# Increases the amount of data points between two points. i.e. inserts 
# 10* the number of points in an array ([0,1] --> [0,0.1,0.2,...,1]
# Used to colour lines in a scatter plot using a cmap. 
#
## truncate_colormap()
# Return a subset of a colormap, handy for sequential cmaps.


def map_formatter2(ax,tick_base_x=15.0, tick_base_y=15.0, labelsize=20,top_label=False, bottom_label=True, right_label=False,left_label=True):
    x_ticks = np.arange(-180,180,tick_base_x)
    y_ticks = np.arange(-90,90,tick_base_y)
    #ax.set_global()
    ax.coastlines(resolution='10m')
    ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.75, color='k', linestyle=':')
    ax.set_xticks(x_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(y_ticks, crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    if top_label:
        ax.xaxis.set_tick_params(labeltop='on',labelsize=labelsize)
    else:
        ax.xaxis.set_tick_params(labeltop='off')
    if right_label:
        ax.yaxis.set_tick_params(labelright='on',labelsize=labelsize)
    else:
        ax.yaxis.set_tick_params(labelright='off')
    if left_label:
        ax.yaxis.set_tick_params(labelleft='on',labelsize=labelsize)
    else:
        ax.yaxis.set_tick_params(labelleft='off')
    if bottom_label:
        ax.xaxis.set_tick_params(labelbottom='on',labelsize=labelsize)
    else:
        ax.xaxis.set_tick_params(labelbottom='off')
    
def map_formatter(ax,tick_base_x=15.0, tick_base_y=15.0, labelsize=20,top_label=False, bottom_label=True, right_label=False,left_label=True,central_longitude=0.0,res='10m'):
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.75, color='k', linestyle=':')
    gl.xlabels_top = top_label
    gl.ylabels_right = right_label
    gl.ylabels_left = left_label
    gl.xlabels_bottom = bottom_label        
    gl.xlocator = mticker.MultipleLocator(base=tick_base_x)
    gl.ylocator = mticker.MultipleLocator(base=tick_base_y)
    ax.coastlines(resolution=res, color='k', linewidth=1)
    gl.xlabel_style= {'size':labelsize}
    gl.ylabel_style= {'size':labelsize}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
def normalise_cmap(vmin, vmax, midpoint, dl, **kwargs):
    if 'cmap' in kwargs:
        act_cmap = kwargs['cmap']
    else:
        act_cmap = cmocean.cm.curl
    levels=np.arange(vmin,vmax,dl)
    midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
    vals = np.interp(midp, [vmin, midpoint, vmax], [0, 0.5, 1])
    colors = act_cmap(vals)
    norm_cmap, norm = from_levels_and_colors(levels, colors)
    return norm_cmap, norm, levels

class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)

def highResPoints(x,y,factor=10):
    '''
    Take points listed in two vectors and return them at a higher
    resultion. Create at least factor*len(x) new points that include the
    original points and those spaced in between.

    Returns new x and y arrays as a tuple (x,y).
    '''
    # r is the distance spanned between pairs of points
    r = [0]
    for i in range(1,len(x)):
        dx = x[i]-x[i-1]
        dy = y[i]-y[i-1]
        r.append(np.sqrt(dx*dx+dy*dy))
    r = np.array(r)

    # rtot is a cumulative sum of r, it's used to save time
    rtot = []
    for i in range(len(r)):
        rtot.append(r[0:i].sum())
    rtot.append(r.sum())

    dr = rtot[-1]/(NPOINTS*RESFACT-1)
    xmod=[x[0]]
    ymod=[y[0]]
    rPos = 0 # current point on walk along data
    rcount = 1 
    while rPos < r.sum():
        x1,x2 = x[rcount-1],x[rcount]
        y1,y2 = y[rcount-1],y[rcount]
        dpos = rPos-rtot[rcount] 
        theta = np.arctan2((x2-x1),(y2-y1))
        rx = np.sin(theta)*dpos+x1
        ry = np.cos(theta)*dpos+y1
        xmod.append(rx)
        ymod.append(ry)
        rPos+=dr
        while rPos > rtot[rcount+1]:
            rPos = rtot[rcount+1]
            rcount+=1
            if rcount>rtot[-1]:
                break
    return xmod,ymod
    
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def extend_cmap(cmap):
    '''
    Ensures the max and min cmap colours are what they're supposed to be
    i.e. sets the cmap up to use "extend='both'" in a contourf plot.
    '''
    min_cmap = cmap(0); max_cmap = cmap(1e6)
    cmap.set_under(min_cmap); cmap.set_over(max_cmap)
    return cmap 

def find_subplot_dims(n):
    n_cols = np.ceil(n**0.5).astype(int)
    n_rows = np.ceil(1.*n / n_cols).astype(int)
    return n_rows,n_cols 
    
    
def fractional_polar_axes(f, thlim=(0, 180), rlim=(0, 1), step=(30, 0.2),
                          thlabel='theta', rlabel='r', ticklabels=True):
    from matplotlib.transforms import Affine2D
    from matplotlib.projections import PolarAxes
    from mpl_toolkits.axisartist import angle_helper
    from mpl_toolkits.axisartist.grid_finder import MaxNLocator
    from mpl_toolkits.axisartist.floating_axes import GridHelperCurveLinear, FloatingSubplot
    """Return polar axes that adhere to desired theta (in deg) and r limits. steps for theta
    and r are really just hints for the locators. Using negative values for rlim causes
    problems for GridHelperCurveLinear for some reason"""
    th0, th1 = thlim # deg
    r0, r1 = rlim
    thstep, rstep = step

    # scale degrees to radians:
    tr_scale = Affine2D().scale(np.pi/180., 1.)
    tr = tr_scale + PolarAxes.PolarTransform()
    theta_grid_locator = angle_helper.LocatorDMS((th1-th0) // thstep)
    r_grid_locator = MaxNLocator((r1-r0) // rstep)
    theta_tick_formatter = angle_helper.FormatterDMS()
    grid_helper = GridHelperCurveLinear(tr,
                                        extremes=(th0, th1, r0, r1),
                                        grid_locator1=theta_grid_locator,
                                        grid_locator2=r_grid_locator,
                                        tick_formatter1=theta_tick_formatter,
                                        tick_formatter2=None)

    a = FloatingSubplot(f, 111, grid_helper=grid_helper)
    f.add_subplot(a)

    # adjust x axis (theta):
    a.axis["bottom"].set_visible(False)
    a.axis["top"].set_axis_direction("bottom") # tick direction
    a.axis["top"].toggle(ticklabels=ticklabels, label=bool(thlabel))
    a.axis["top"].major_ticklabels.set_axis_direction("top")
    a.axis["top"].label.set_axis_direction("top")

    # adjust y axis (r):
    a.axis["left"].set_axis_direction("bottom") # tick direction
    a.axis["right"].set_axis_direction("top") # tick direction
    a.axis["left"].toggle(ticklabels=ticklabels, label=bool(rlabel))

    # add labels:
    a.axis["top"].label.set_text(thlabel)
    a.axis["left"].label.set_text(rlabel)

    # create a parasite axes whose transData is theta, r:
    auxa = a.get_aux_axes(tr)
    # make aux_ax to have a clip path as in a?:
    auxa.patch = a.patch 
    # this has a side effect that the patch is drawn twice, and possibly over some other
    # artists. So, we decrease the zorder a bit to prevent this:
    a.patch.zorder = -2

    # add sector lines for both dimensions:
    thticks = grid_helper.grid_info['lon_info'][0]
    rticks = grid_helper.grid_info['lat_info'][0]
    for th in thticks[1:-1]: # all but the first and last
        auxa.plot([th, th], [r0, r1], '--', c='grey', zorder=-1)
    for ri, r in enumerate(rticks):
        # plot first r line as axes border in solid black only if it isn't at r=0
        if ri == 0 and r != 0:
            ls, lw, color = 'solid', 2, 'black'
        else:
            ls, lw, color = 'dashed', 1, 'grey'
        # From http://stackoverflow.com/a/19828753/2020363
        auxa.add_artist(plt.Circle([0, 0], radius=r, ls=ls, lw=lw, color=color, fill=False,
                        transform=auxa.transData._b, zorder=-1))
    return auxa


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def reverse_colourmap(cmap, name = 'my_cmap_r'):
    """
    In: 
    cmap, name 
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
    """        
    reverse = []
    k = []   

    for key in cmap._segmentdata:    
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:                    
            data.append((1-t[0],t[2],t[1]))            
        reverse.append(sorted(data))    

    LinearL = dict(zip(k,reverse))
    my_cmap_r = matplotlib.colors.LinearSegmentedColormap(name, LinearL) 
    return my_cmap_r
