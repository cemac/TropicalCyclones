import iris
import iris.analysis.calculus
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from plotting_tools import *
from storm_tools import *
from cube_tools import *
import matplotlib.cm as mpl_cm
import iris.plot as iplt

def load_ens_members(fpath, dd, hr, yr=2014, mth=12):
    md  = 1203; TT = 12
    NN = (dd - (md-1200))*24 + (hr - TT)
    print NN
    u = []
    v = []
    # fpath should be everything but ens_member.pp
    ensemble_members=np.arange(12)
    stash_code1 = 'm01s15i201'
    stash_code2 = 'm01s15i202'
    #plev = 850 #Set to 0 for single level data
    min_val = 0
    max_val = 0
    u_list = []; v_list = []
    data_constraint1 = iris.AttributeConstraint(STASH=stash_code1)
    data_constraint2 = iris.AttributeConstraint(STASH=stash_code2)
    cenlats = []; cenlons = []
#    p_constraint = iris.Constraint(pressure=plev)
    time_constraint = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr))
    for em in ensemble_members:
        [cenlat,cenlon] = find_centre_ens_glm(md,TT,em,NN)
        cenlats.append(cenlat)
        cenlons.append(cenlon)
        minlat = cenlat-10;maxlat = cenlat+10; minlon = cenlon-10; maxlon = cenlon+10;
        b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
        df = fpath + '{0:02d}.pp'.format(em)
        lat_sample = [('latitude',cenlat)]
        lon_sample = [('longitude',cenlon)]
        with iris.FUTURE.context(cell_datetime_objects=True):
            u = iris.load_cube(df,data_constraint1).extract(time_constraint&b_constraint)
            v = iris.load_cube(df,data_constraint2).extract(time_constraint&b_constraint)
        u_lons = u.interpolate(lat_sample,iris.analysis.Linear())
        v_lons = v.interpolate(lat_sample,iris.analysis.Linear())
        #v_lats = data.interpolate(lon_sample,iris.analysis.Linear())
        u_list.append(u_lons)
        v_list.append(v_lons)
    return u_list,v_list, cenlats, cenlons
    
def plot_ws_slc(u_list,v_list,cenlat,cenlon,outfile):
    fig = plt.figure(figsize=(20,10), dpi=100)
    plt.gcf().subplots_adjust(hspace=0.05, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925)
    v_max = 0; v_min = 0
    ens_members = np.arange(12)
    
    ws_max = 30
    ws_min = -30
    dl = 5
    brewer_cmap = mpl_cm.get_cmap('brewer_RdBu_11')
    
    maxcmap = brewer_cmap(1e7); mincmap = brewer_cmap(0)
    
    cmap, norm, ws_levels = normalise_cmap(ws_min,ws_max+1,0,dl,cmap=brewer_cmap)
    
    
    for em in ens_members:
        ax = plt.subplot(3,4,em+1)
        v = v_list[em]
        x = v.coord('longitude').points - cenlon[em]
        z = v.coord('pressure').points
        X,Z = np.meshgrid(x,z)
        cs = ax.contourf(X,Z,v.data,levels=ws_levels,norm=norm,cmap=cmap,extend='both')
        cs.cmap.set_under(mincmap); cs.cmap.set_over(maxcmap)
        ax.invert_yaxis()
        if em % 4 !=0:
            ax.set_yticklabels([])
        if em < 8:
            ax.set_xticklabels([])
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.05])
    colorbar = plt.colorbar(cs, cbar_ax, orientation='horizontal')
    
    #plt.show()
    plt.savefig(outfile)
    plt.close()

def main():
    monthday = [1203]; times = [12]
    vt_day = [3,4,5,7,8]; vt_hr = [0, 6, 12, 18]
    #### Load cubes:
    for md in monthday:
        for TT in times:
            dt = 0
            for dd in vt_day:
                for hr in vt_hr:
                    if dd==(md-1200) and hr < TT:
                        continue
                    outfile = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/plots/4p4/{0:04d}_{1:02d}Z/ensemble_plots/vortex_cs_{2:02d}_{3:02d}Z.png'.format(md,TT,dd,hr)
                    fpath = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/4p4/{0:04d}_{1:02d}Z/wind_plev/uvplev_4p4_{0:04d}_{1:02d}Z_em'.format(md,TT)
                    u_list, v_list, cenlat, cenlon = load_ens_members(fpath,dd,hr)
                    plot_ws_slc(u_list,v_list,cenlat, cenlon,outfile)
                    dt = dt + 6
                    
if __name__=='__main__':
    main()
