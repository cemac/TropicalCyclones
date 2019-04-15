########################################################################
# Stamp plots of wind on pressure levels. Filled contour and quiver.   #
#                                          John Ashcroft, October 2017 #
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

def load_ens_members(stash,fpath, dd, hr,plev=0, yr=2014, mth=12):
    # fpath should be everything but ens_member.pp
    ensemble_members=np.arange(12)
    minlat = -10;maxlat = 40; minlon = 80; maxlon = 160;
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    ens_members = np.arange(12)
    cube_list = [];
    data_constraint = iris.AttributeConstraint(STASH=stash)
    p_constraint = iris.Constraint(pressure=plev)
    time_constraint = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr))
        
    for em in ensemble_members:
        df = fpath + '{0:02d}.pp'.format(em)
        if plev == 0:
            with iris.FUTURE.context(cell_datetime_objects=True):
                data = iris.load_cube(df,data_constraint).extract(time_constraint&b_constraint)
        else:
            with iris.FUTURE.context(cell_datetime_objects=True):
                data = iris.load_cube(df,data_constraint).extract(p_constraint&time_constraint&b_constraint)
        cube_list.append(data); 
    return cube_list

def myround(x, base=5):
    return int(base * round(float(x)/base))
    
def load_era_data(era_name,cube,dd,hr,plev,yr=2014,mth=12):
    #plev=850
    TT = myround(hr,base=6)
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    minlat = np.amin(lats); maxlat = np.amax(lats); minlon = np.amin(lons); maxlon = np.amax(lons)
    df = '/nfs/a137/earceb/ERA-interim/{0:04d}/as/ggas{0:04d}{1:02d}{2:02d}{3:02d}00.nc'.format(yr, mth, dd,TT)
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    if plev == 0:
        data = iris.load_cube(df,era_name).extract(b_constraint)[0][:][:]
    else:
        p_constraint = iris.Constraint(p=plev)
        data = iris.load_cube(df,era_name).extract(b_constraint&p_constraint)[0][:][:]
    return data
    
def plot_data(data_list,avg_data,outfile,md,TT,pos_id,dtick=5):
    #brewer_cmap = mpl_cm.get_cmap('brewer_OrRd_09')
    ensemble_members=np.arange(12)
    cmap = mpl_cm.get_cmap('brewer_Accent_08')
    
    for em in ensemble_members:
        dataname = '/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver_glm_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(md,TT,em)
        track_data = np.load(dataname)
        if em == 0:
            lats = track_data['lats']
            lons = track_data['lons']
        else:
            lats = np.vstack([lats,track_data['lats']])
            lons = np.vstack([lons,track_data['lons']])
    
    
    for em in ensemble_members:
        cube = data_list[em];
        cube = cube - avg_data
        min_v = np.amin(cube.data)
        max_v = np.amax(cube.data)
        if em == 0: 
            min_val = min_v; max_val = max_v
        else:
            if min_v < min_val:
                min_val = min_v
            if max_v > max_val:
                max_val = max_v
    
    #min_val = np.floor_divide(min_val,1) * 1
    #max_val = np.floor_divide(max_val,1) * 1
    #dl = (max_val - min_val)/20.
    brewer_cmap = mpl_cm.get_cmap('brewer_Blues_09')
    brewer_cmap = truncate_colormap(brewer_cmap,0.0,0.7)
    #maxcmap = brewer_cmap(1e7); mincmap = brewer_cmap(0)
    #cmap, norm, cube_levs = normalise_cmap(min_val,max_val,0,dl,cmap=brewer_cmap)
    
    cloud_levs = np.arange(0,1.00001,0.05)
    
    
    
    plt.figure(figsize=(20,10), dpi=100)
    plt.gcf().subplots_adjust(hspace=0.05, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925)
    for em in ensemble_members:
        cube = data_list[em];
        #cube = cube - avg_data
        plt.subplot(4,4,em+1,projection=ccrs.PlateCarree())
        cube_plt = iplt.contourf(cube,levels=cloud_levs,cmap=brewer_cmap)
#        ws_plt.cmap.set_under('b'); ws_plt.cmap.set_over('r')
        ax = plt.gca()
        if em % 4 != 0:
            left_label=False
        else:
            left_label=True
        if em < 9 or em == 11:
            bottom_label=False
        else:
            bottom_label=True
        for i in ensemble_members:
            if i!=em:
                ax.plot(lons[i,:],lats[i,:],color='lightgrey',alpha=0.4)
            else:
                mem = i
        ax.plot(lons[mem,:],lats[mem,:],color='m')
        ax.plot(lons[mem,pos_id],lats[mem,pos_id],'r*')
        map_formatter(ax,bottom_label=bottom_label, left_label=left_label, labelsize=10,tick_base_x=10,tick_base_y=10)
        ax.annotate('{0:02d}'.format(em), xy=(0.99,0.01), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white',fontsize=20)
        
    
    #plt.subplot(4,4,13)
    #T_plt = iplt.contourf(avg_T,levels=T_levs,cmap='plasma')
    #ws_plt.cmap.set_under('b'); ws_plt.cmap.set_over('r')
    #ax = plt.gca()
    #map_formatter(ax,labelsize=12,tick_base_x=10,tick_base_y=10)
    #ax.annotate('avg', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
        
    ### ERA data
    #plt.subplot(4,4,16)
    #T_plt = iplt.contourf(era_T,levels=T_levs,cmap='plasma')
    #    ws_plt.cmap.set_under('b'); ws_plt.cmap.set_over('r')
    #ax = plt.gca()
    #map_formatter(ax,labelsize=12,tick_base_x=10,tick_base_y=10,left_label=False)
    #ax.annotate('ERA', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
    
    colorbar_axes = plt.gcf().add_axes([0.35, 0.05, 0.3, 0.05])
    colorbar = plt.colorbar(cube_plt, colorbar_axes, orientation='horizontal')
    plt.suptitle('Low clouds')    
    #time_unit, validity_time, data_time, lead_time = get_time(avg_data)
    #ax.annotate('Simulation start time: {:%Y/%m/%d T%H:%MZ}, Lead time: +{:02d}Z'.format(data_time, int(lead_time)),xy=(0.35,0.2),horizontalalignment='left', verticalalignment='top',xycoords='figure fraction')
    #iplt.show()
    plt.savefig(outfile)
    plt.close()
    
    
def main():
    monthday = [1203]; times = [12]
    vt_day = [3,4,5,6,7]; vt_hr = [0, 6, 12, 18]
    plevs = [0]; # 0 for single layer data
    stash = 'm01s09i203'; # MED cloud (Low = 203, Med = 204, High = 205)
    #### Load cubes:
    for md in monthday:
        for TT in times:
            pos_id = 0
            for dd in vt_day:
                for hr in vt_hr:
                    if dd==(md-1200) and hr < TT:
                        continue
                    for plev in plevs:
                        outfile = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/plots/4p4/{0:04d}_{1:02d}Z/ensemble_plots/lowclouds_{2:02d}_{3:02d}Z_{4:03d}hPa.png'.format(md,TT,dd,hr,plev)
                        fpath = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/glm/{0:04d}_{1:02d}Z/clouds/clouds_glm_{0:04d}_{1:02d}Z_em'.format(md,TT)
                        data_list = load_ens_members(stash,fpath,dd,hr,plev=plev)
                        avg_data = avg_cubes(data_list); 
                        plot_data(data_list,avg_data,outfile,md,TT,pos_id)
                        pos_id = pos_id + 1
                    
                    


        
    
if __name__=='__main__':
    main()
