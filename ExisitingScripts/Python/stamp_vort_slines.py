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

from windspharm.iris import VectorWind


def load_ens_members(md,TT,dd,hr,mth=12,yr=2014):
    vort_plev = 200
    sline_plev = 200
    minlat = -10;maxlat = 30; minlon = 115; maxlon = 180;
    ens_members = np.arange(12)
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    fpath = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/glm/{0:04d}_{1:02d}Z/'.format(md,TT)
    wind_ext = 'wind_plev/uvplev_glm_{0:04d}_{1:02d}Z_'.format(md,TT)
    
    u_stash = 'm01s15i201'
    v_stash = 'm01s15i202'
        
    uvort_list = []; vvort_list = []; usline_list = []; vsline_list = []
    u_constraint = iris.AttributeConstraint(STASH=u_stash)
    v_constraint = iris.AttributeConstraint(STASH=v_stash)
    vortp_constraint = iris.Constraint(pressure=vort_plev)
    slinep_constraint = iris.Constraint(pressure=sline_plev)
    time_constraint = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr))
    
    for em in ens_members:
        wind_df = fpath + wind_ext + 'em{0:02d}.pp'.format(em)
        with iris.FUTURE.context(cell_datetime_objects=True):
            uvort = iris.load_cube(wind_df,u_constraint).extract(vortp_constraint&time_constraint)
            vvort = iris.load_cube(wind_df,v_constraint).extract(vortp_constraint&time_constraint)
            usline = iris.load_cube(wind_df,u_constraint).extract(slinep_constraint&time_constraint)
            vsline = iris.load_cube(wind_df,v_constraint).extract(slinep_constraint&time_constraint)
        uvort_list.append(uvort); vvort_list.append(vvort); vsline_list.append(vsline); usline_list.append(usline)
    return uvort_list, vvort_list, usline_list, vsline_list 
        
    
def plot_data(vort,sline,avg_vort,avg_sline,outfile,pos_id):
    vrtmax = 21e-5
    vrtmin = -18e-5
    dl = 3e-5
    md = 1203; TT = 12
    brewer_cmap = mpl_cm.get_cmap('brewer_RdBu_11')
    maxcmap = brewer_cmap(1e7); mincmap = brewer_cmap(0)
    cmap, norm, vortlevels = normalise_cmap(vrtmin,vrtmax,0,dl,cmap=brewer_cmap)
    psi_levs = np.arange(-1e8,1e8,6e6)
    
    ens_members=np.arange(12)
    for em in ens_members:
        dataname = '/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver_glm_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(1203,12,em)
        track_data = np.load(dataname)
        if em == 0:
            lats = track_data['lats']
            lons = track_data['lons']
        else:
            lats = np.vstack([lats,track_data['lats']])
            lons = np.vstack([lons,track_data['lons']])
    
    
    fig=plt.figure(figsize=(20,10), dpi=100)
    plt.gcf().subplots_adjust(hspace=0.05, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925)
    for em in ens_members:
        vrt = vort[em];
        psi = sline[em]
        
        plt.subplot(3,4,em+1,projection=ccrs.PlateCarree())
        cube_plt = iplt.contourf(vrt,cmap=cmap,levels=vortlevels,norm=norm,extend='both')
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        psi_plt = iplt.contour(psi,levels=psi_levs,colors='darkcyan',linewidths=2)
        cube_plt.cmap.set_under(mincmap); cube_plt.cmap.set_over(maxcmap)
        ax = plt.gca()
        if em % 4 != 0:
            left_label=False
        else:
            left_label=True
        if em < 8:
            bottom_label=False
        else:
            bottom_label=True
        for i in ens_members:
            if i!=em:
                ax.plot(lons[i,:],lats[i,:],color='lightgrey',alpha=0.4)
            else:
                mem = i
        ax.plot(lons[mem,:],lats[mem,:],color='b')
        ax.plot(lons[mem,pos_id],lats[mem,pos_id],'r*')
        map_formatter(ax,bottom_label=bottom_label, left_label=left_label, labelsize=10,tick_base_x=10,tick_base_y=10)
        ax.annotate('{0:02d}'.format(em), xy=(0.99,0.01), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white',fontsize=20)
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.05])
    colorbar = plt.colorbar(cube_plt, cbar_ax, orientation='horizontal',format=OOMFormatter(-5, mathText=False))
   
    print outfile
    
    #plt.show()
    plt.savefig(outfile)
    plt.close()
    

def main():
    monthday = [1203]; times = [12]
    vt_day = [3,4,5,6,7]; vt_hr = [0, 6, 12, 18]
    minlat = -10;maxlat = 30; minlon = 115; maxlon = 180;
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    psi_list = []
    #### Load cubes:
    for md in monthday:
        for TT in times:
            pos_id = 0
            for dd in vt_day:
                for hr in vt_hr:
                    vort = []; sline = []
                    if dd==(md-1200) and hr < TT:
                        continue
                    outfile = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/plots/glm/1203_12Z/story_time/stamp_vrt200_sline200_{0:02d}_{1:02d}Z.png'.format(dd,hr)
                    uvort_list, vvort_list, usline_list, vsline_list = load_ens_members(md,TT,dd,hr)
                    avg_uvort = avg_cubes(uvort_list); avg_vvort = avg_cubes(vvort_list); 
                    avg_usline = avg_cubes(usline_list); avg_vsline = avg_cubes(vsline_list); 
                    wvort = VectorWind(avg_uvort,avg_vvort); wsline = VectorWind(avg_usline,avg_vsline)
                    avg_sline = wsline.streamfunction()
                    avg_vort = wvort.vorticity()
                    for u, v in zip(uvort_list, vvort_list):
                        w = VectorWind(u,v)
                        vrt = w.vorticity()
                        vrt = vrt.extract(b_constraint)
                        vort.append(vrt)
                    for u,v in zip(usline_list,vsline_list):
                        w = VectorWind(u,v)
                        sfunc = w.streamfunction()
                        sfunc = sfunc.extract(b_constraint)
                        sline.append(sfunc)
                        

                    avg_vort = avg_vort.extract(b_constraint)
                    avg_sline = avg_sline.extract(b_constraint)
                    plot_data(vort,sline,avg_vort,avg_sline,outfile,pos_id)
                    pos_id = pos_id + 1

if __name__=='__main__':
    main()
