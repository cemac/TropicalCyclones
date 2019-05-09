########################################################################
# SEnsemble average of winds and geopotential.                         #
#                                         John Ashcroft, November 2017 #
########################################################################


import iris
import iris.analysis.calculus
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from plotting_tools import *
from cube_tools import *
import matplotlib.cm as mpl_cm
import iris.plot as iplt


def load_ens_members(md,TT,dd,hr,mth=12,yr=2014):
    geop_plev = 500
    wind_plev = 500
    minlat = -10;maxlat = 30; minlon = 115; maxlon = 180;
    ens_members = np.arange(12)
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    fpath = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/glm/{0:04d}_{1:02d}Z/'.format(md,TT)
    geop_ext = 'tctracker/tcver_glm_{0:04d}_{1:02d}Z_'.format(md,TT)
    wind_ext = 'tctracker/tcver_glm_{0:04d}_{1:02d}Z_'.format(md,TT)
    
#    geop_stash = 'm01s16i202'
#    u_stash = 'm01s15i201'
#    v_stash = 'm01s15i202'

    geop_stash = 'm01s16i222'
    u_stash = 'm01s03i225'
    v_stash = 'm01s03i226'
        
    u_list = []; v_list = []; geop_list = []
    u_constraint = iris.AttributeConstraint(STASH=u_stash)
    v_constraint = iris.AttributeConstraint(STASH=v_stash)
    geop_constraint = iris.AttributeConstraint(STASH=geop_stash)
    geopp_constraint = iris.Constraint(pressure=geop_plev)
    windp_constraint = iris.Constraint(pressure=wind_plev)
    time_constraint = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr+1))
    
    for em in ens_members:
        geop_df = fpath + geop_ext + 'em{0:02d}.pp'.format(em)
        wind_df = fpath + wind_ext + 'em{0:02d}.pp'.format(em)
        with iris.FUTURE.context(cell_datetime_objects=True):
            u = iris.load_cube(wind_df,u_constraint).extract(time_constraint&b_constraint)
            v = iris.load_cube(wind_df,v_constraint).extract(time_constraint&b_constraint)
            geop = iris.load_cube(geop_df,geop_constraint).extract(time_constraint&b_constraint)
        u_list.append(u); v_list.append(v); geop_list.append(geop)
    return u_list, v_list, geop_list
        
    
def plot_data(u,v,geop,outfile):
    step = 10; scale = 300
    max_ws = 80
    ws = (u**2 + v**2)**0.5
    ws.rename('windspeed')
    
    mingeop = np.floor_divide(np.amin(geop.data),40) * 40
    maxgeop = (np.floor_divide(np.amax(geop.data),40) + 2) * 40
    #geoplevels = np.arange(mingeop,maxgeop,40)
    
    geop.data = geop.data/100.
    minslp = np.floor_divide(np.amin(geop.data),4)*4
    maxslp = (np.floor_divide(np.amax(geop.data),4) +1)*4
    geoplevels = np.arange(minslp,maxslp,4)
    
    cmap = mpl_cm.get_cmap('plasma')
    max_wscmap = cmap(1e7)
    ws_levs = np.arange(0,max_ws,3)
    
    x = u.coord('longitude').points
    y = u.coord('latitude').points
    
    
    fig = plt.figure(1,figsize=(20,10),)
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    avg_ws_plt = iplt.contourf(ws,levels=ws_levs,cmap=cmap,extend='max')
    avg_ws_plt.cmap.set_over(max_wscmap)
    plt.quiver(x[::step],y[::step],u.data[::step,::step],v.data[::step,::step],pivot='middle',transform=ccrs.PlateCarree(),color='w',scale=scale)
    geop_plt = iplt.contour(geop,levels=geoplevels,colors='k',linewidths=2)
    plt.clabel(geop_plt,fmt = '%1.0f')
    map_formatter(ax,bottom_label=True, left_label=True, labelsize=15,tick_base_x=15,tick_base_y=10)
    cbar = plt.colorbar(avg_ws_plt,orientation='vertical',fraction=0.046,pad=0.09)
    
#    plt.show()
    plt.savefig(outfile)
    plt.close()
    

def main():
    monthday = [1203]; times = [12]
    vt_day = [3,4,5,6,7]; vt_hr = [0, 6, 12, 18]
    #### Load cubes:
    for md in monthday:
        for TT in times:
            for dd in vt_day:
                for hr in vt_hr:
                    if dd==(md-1200) and hr < TT:
                        continue
                    outfile = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/plots/glm/1203_12Z/story_time/geop_wind_b_{0:02d}_{1:02d}Z.png'.format(dd,hr)
                    u_list, v_list, geop_list = load_ens_members(md,TT,dd,hr)
                    avg_u = avg_cubes(u_list); avg_v = avg_cubes(v_list); avg_geop = avg_cubes(geop_list)
                    plot_data(avg_u,avg_v,avg_geop,outfile)

if __name__=='__main__':
    main()
