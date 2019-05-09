########################################################################
# Script to produce stamp plots of vorticity and sea level pressure.   #
#                                              John Ashcroft, Oct 2017 #
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
from vorticity import calc_vrt_spherical

def load_ens_members(fpath,slp_fpath, dd, hr, yr=2014, mth=12):
    # Load each of the 12 ensemble members. Return a cube list of vorticity
    # and slp.
    # fpath should be everything but ens_member.pp
    
    # Define the domain for plotting:- 
    minlat = 0; maxlat = 45; minlon = 105; maxlon = 180
    
    # 12 Ensemble members, define the stash codes and pressure levels for required data
    # Make iris constraints for loading cubes
    ensemble_members=np.arange(12)
    stash_code1 = 'm01s15i201'
    stash_code2 = 'm01s15i202'
    slp_stash = 'm01s16i222'
    p_vrt = 850 #Set to 0 for single level data
    p_shear_high = 200; p_shear_low = 850
    ens_members = np.arange(12)
    min_val = 0
    max_val = 0
    vrt_list = []; ushear_list = []; vshear_list = []; slp_list = []
    data_constraint1 = iris.AttributeConstraint(STASH=stash_code1)
    data_constraint2 = iris.AttributeConstraint(STASH=stash_code2)
    slp_constraint = iris.AttributeConstraint(STASH=slp_stash)
    pvrt_constraint = iris.Constraint(pressure=p_vrt)
    pshearh_constraint = iris.Constraint(pressure=p_shear_high)
    pshearl_constraint = iris.Constraint(pressure=p_shear_low)
    time_constraint = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr))
    time_constraint2 = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr+1))
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    
    ## Load all of the cubes.        
    for em in ensemble_members:
        df = fpath + '{0:02d}.pp'.format(em)
        slp_df = slp_fpath + '{0:02d}.pp'.format(em)
        with iris.FUTURE.context(cell_datetime_objects=True):
            uvrt = iris.load_cube(df,data_constraint1).extract(pvrt_constraint&time_constraint&b_constraint)
            vvrt = iris.load_cube(df,data_constraint2).extract(pvrt_constraint&time_constraint&b_constraint)
            uhigh = iris.load_cube(df,data_constraint1).extract(pshearh_constraint&time_constraint&b_constraint)
            vhigh = iris.load_cube(df,data_constraint2).extract(pshearh_constraint&time_constraint&b_constraint)
            ulow = iris.load_cube(df,data_constraint1).extract(pshearl_constraint&time_constraint&b_constraint)
            vlow = iris.load_cube(df,data_constraint2).extract(pshearl_constraint&time_constraint&b_constraint)
            slp = iris.load_cube(slp_df,slp_constraint).extract(time_constraint2&b_constraint)
            vrt = iris.analysis.calculus.curl(uvrt,vvrt)[2]
            slp.data = slp.data/100.
            ushear = uhigh-ulow; vshear=vhigh-vlow
        vrt_list.append(vrt); slp_list.append(slp); ushear_list.append(ushear); vshear_list.append(vshear)
        #slp_list.append(slp)
    return vrt_list,ushear_list,vshear_list,slp_list

def myround(x, base=5):
    return int(base * round(float(x)/base))
    
def load_era_data(cube,dd,hr,yr=2014,mth=12):
    p_vrt = 850 
    p_shear_high = 200; p_shear_low = 850
    TT = myround(hr,base=6)
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    minlat = np.amin(lats); maxlat = np.amax(lats); minlon = np.amin(lons); maxlon = np.amax(lons)
    df = '/nfs/a137/earceb/ERA-interim/{0:04d}/ap/ggap{0:04d}{1:02d}{2:02d}{3:02d}00.nc'.format(yr, mth, dd,TT)
    pvrt_constraint = iris.Constraint(p=p_vrt)
    phigh_constraint = iris.Constraint(p=p_shear_high)
    plow_constraint = iris.Constraint(p=p_shear_low)
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    
    u = iris.load_cube(df,'eastward_wind')[0][:][:][:]
    v = iris.load_cube(df,'northward_wind')[0][:][:][:]
    uhigh = u.extract(phigh_constraint&b_constraint)
    ulow = u.extract(plow_constraint&b_constraint)
    vhigh = v.extract(phigh_constraint&b_constraint)
    vlow = v.extract(plow_constraint&b_constraint)
    era_ushear = uhigh - ulow
    era_vshear = vhigh - vlow
    print iris.load(df)
    #exit()
    era_vrt = iris.load(df)[6].extract(pvrt_constraint&b_constraint)[0][:][:]  # Vorticity
    
    df = '/nfs/a137/earceb/ERA-interim/{0:04d}/as/ggas{0:04d}{1:02d}{2:02d}{3:02d}00.nc'.format(yr, mth, dd,TT)
    era_slp = iris.load_cube(df,'air_pressure_at_sea_level').extract(b_constraint)[0][0][:][:]
    era_slp.data = era_slp.data / 100.
    return era_vrt, era_ushear,era_vshear,era_slp
    
    
def round_down(num, divisor):
    return num - (num%divisor)

def round_up(num,divisor):
    return num + divisor - (num%divisor)

    
def plot_data(vrt_list,ushear_list,vshear_list,slp_list,era_vrt,era_ushear,era_vshear,era_slp,outfile,dtick=5):
    #brewer_cmap = mpl_cm.get_cmap('brewer_OrRd_09')
    scale = 500.; step=10
    vrtmax = 18e-5
    vrtmin = -15e-5
    dl = 3e-5
    md = 1203; TT = 12
    brewer_cmap = mpl_cm.get_cmap('brewer_RdBu_11')
    
    maxcmap = brewer_cmap(1e7); mincmap = brewer_cmap(0)
    shear_levs = np.arange(0,50,5)
    #ens_members = np.arange(12)
    #for i in ens_members:
    #    if i == 0:
    #        vrtmin = np.amin(vrt_list[i].data)
    #        vrtmax = np.amax(vrt_list[i].data)
    #    else:
    #        max_val = np.amax(vrt_list[i].data)
    #        min_val = np.amin(vrt_list[i].data)
    #        if max_val > vrtmax:
    #            vrtmax = max_val
    #        if min_val < vrtmin:
    #            vrtmin = min_val
    
    cmap, norm, vortlevels = normalise_cmap(vrtmin,vrtmax,0,dl,cmap=brewer_cmap)
    
    avg_vrt = avg_cubes(vrt_list)
    avg_ushear = avg_cubes(ushear_list)
    avg_vshear = avg_cubes(vshear_list)
    avg_slp = avg_cubes(slp_list)
    
    ensemble_members=np.arange(12)
    
    for em in ensemble_members:
        dataname = '/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver_glm_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(md,TT,em)
        track_data = np.load(dataname)
        if em == 0:
            lats = track_data['lats']
            lons = track_data['lons']
        else:
            lats = np.vstack([lats,track_data['lats']])
            lons = np.vstack([lons,track_data['lons']])

    fig = plt.figure(figsize=(20,10), dpi=100)
    plt.gcf().subplots_adjust(hspace=0.05, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925)
    #plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    cmap.set_under(mincmap); cmap.set_over(maxcmap)
        
    for em in ensemble_members:
        vrt = vrt_list[em]; ushear = ushear_list[em]
        slp = slp_list[em]; vshear = vshear_list[em]
        #vrt.data[:][:] = vrt.data[:][:] - avg_vrt.data[:][:]
        maxslp = round_up(np.amax(slp.data),4)
        minslp = round_down(np.amin(slp.data),4)
        slplevels = np.arange(minslp,maxslp,4)
        plt.subplot(3,4,em+1)
        #slp = iplt.contour(slp,levels=slplevels, colors='darkslategrey')
        #plt.clabel(slp)
        vrt = iplt.contourf(vrt,levels=vortlevels,cmap=cmap,norm=norm,extend='both')
        #if em == 0:
        #    vrt = iplt.contourf(vrt,levels=vortlevels,cmap=cmap,norm=norm,extend='both')
        #    vrt.cmap.set_under(mincmap); vrt.cmap.set_over(maxcmap)
        
        x = ushear.coord('longitude').points
        y = ushear.coord('latitude').points
        X,Y = np.meshgrid(x,y)
        ax = plt.gca()
        ax.set_ylim([0,45])
        arrows = ax.quiver(X[::step, ::step],Y[::step,::step],ushear.data[::step,::step],vshear.data[::step,::step], angles='xy', scale=scale, transform=ccrs.PlateCarree())
        if em == 2:
            ax.quiverkey(arrows, 0.90, 0.97, 25, '25 ms$^{-1}$',labelpos='E',
                   coordinates='figure')
        if em % 4 != 0:
            left_label=False
        else:
            left_label=True
        if em < 8:
            bottom_label=False
        else:
            bottom_label=True
        map_formatter(ax,bottom_label=bottom_label, left_label=left_label, labelsize=11,tick_base_x=10,tick_base_y=10)
        ax.annotate('{0:02d}'.format(em), xy=(0.99,0.01), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white',fontsize=20)
        for i in ensemble_members:
            if i!=em:
                ax.plot(lons[i,:],lats[i,:],color='lightgrey',alpha=0.4)
            else:
                mem = i
        ax.plot(lons[mem,:],lats[mem,:],color='b')
    
    #plt.subplot(4,4,13)
    #maxslp = round_up(np.amax(avg_slp.data),4)
    #minslp = round_down(np.amin(avg_slp.data),4)
    #slplevels = np.arange(minslp,maxslp,4)
    #slp = iplt.contour(avg_slp,levels=slplevels,colors='k')
    #vrt = iplt.contourf(avg_vrt,levels=vortlevels,cmap=cmap,norm=norm,extend='both')
    #vrt.cmap.set_under(mincmap); vrt.cmap.set_over(maxcmap)
    #x = avg_ushear.coord('longitude').points
    #y = avg_ushear.coord('latitude').points
    #X,Y = np.meshgrid(x,y)
    #ax = plt.gca()
    #arrows = ax.quiver(X[::step, ::step],Y[::step,::step],avg_ushear.data[::step,::step],avg_vshear.data[::step,::step], angles='xy', scale=scale, transform=ccrs.PlateCarree())
    #map_formatter(ax,labelsize=12,tick_base_x=10,tick_base_y=10)
    #ax.annotate('avg', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white',fontsize=20)
        
    #### ERA data
    #plt.subplot(4,4,16)
    #vrt = iplt.contourf(era_vrt,levels=vortlevels,cmap=cmap,norm=norm,extend='both')
    #vrt.cmap.set_under(mincmap); vrt.cmap.set_over(maxcmap)
    #maxslp = round_up(np.amax(era_slp.data),4) 
    #minslp = round_down(np.amin(era_slp.data),4)
    #slplevels = np.arange(minslp,maxslp,4)
    #slp = iplt.contour(era_slp,levels=slplevels,colors='k')
    #x = era_ushear.coord('longitude').points
    #y = era_ushear.coord('latitude').points
    #X,Y = np.meshgrid(x,y)
    #ax = plt.gca()
    #arrows = ax.quiver(X[::step, ::step],Y[::step,::step],era_ushear.data[::step,::step],era_vshear.data[::step,::step], angles='xy', scale=scale, transform=ccrs.PlateCarree())
    #map_formatter(ax,labelsize=12,tick_base_x=10,tick_base_y=10,left_label=False)
    #ax.annotate('ERA-Interim', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white',fontsize=20)
    
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.05])
    colorbar = plt.colorbar(vrt, cbar_ax, orientation='horizontal',format=OOMFormatter(-5, mathText=False))
    colorbar.ax.tick_params(labelsize=14)

    #plt.suptitle('850-200hPa windshear, 850hPa relative vorticity, sea level pressure')    
    time_unit, validity_time, data_time, lead_time = get_time(avg_vrt)
    #ax.annotate('Simulation start time: {:%Y/%m/%d T%H:%MZ}, Lead time: +{:02d}Z'.format(data_time, int(lead_time)),xy=(0.35,0.2),horizontalalignment='left', verticalalignment='top',xycoords='figure fraction')
    #iplt.show()
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
                    print 'dd={0:02d}, hr={1:02d}'.format(dd,hr)
                    outfile = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/plots/glm/{0:04d}_{1:02d}Z/ensemble_plots/slp_vorticityanom_shear_{2:02d}_{3:02d}Z.png'.format(md,TT,dd,hr)
                    fpath = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/glm/{0:04d}_{1:02d}Z/wind_plev/uvplev_glm_{0:04d}_{1:02d}Z_em'.format(md,TT)
                    slp_fpath = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/glm/{0:04d}_{1:02d}Z/tctracker/tcver_glm_{0:04d}_{1:02d}Z_em'.format(md,TT)
                    vrt_list,ushear_list,vshear_list,slp_list = load_ens_members(fpath,slp_fpath,dd,hr)
                    
                    era_vrt, era_ushear,era_vshear,era_slp = load_era_data(vrt_list[0],dd,hr)
                    plot_data(vrt_list,ushear_list,vshear_list,slp_list,era_vrt,era_ushear,era_vshear,era_slp,outfile)
                    
                    


        
    
if __name__=='__main__':
    main()
