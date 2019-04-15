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
import matplotlib.image as mpimg

def load_ens_members(fpath, dd, hr,plev, yr=2014, mth=12):
    # fpath should be everything but ens_member.pp
    minlat = -10;maxlat = 30; minlon = 105; maxlon = 180;
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    ensemble_members=np.arange(12)
    stash_code1 = 'm01s15i201'
    stash_code2 = 'm01s15i202'
    ens_members = np.arange(12)
    min_val = 0
    max_val = 0
    psi_list = []
    data_constraint1 = iris.AttributeConstraint(STASH=stash_code1)
    data_constraint2 = iris.AttributeConstraint(STASH=stash_code2)
    p_constraint = iris.Constraint(pressure=plev)
    time_constraint = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr))
        
    for em in ensemble_members:
        df = fpath + '{0:02d}.pp'.format(em)
        with iris.FUTURE.context(cell_datetime_objects=True):
            u = iris.load_cube(df,data_constraint1).extract(time_constraint)
            v = iris.load_cube(df,data_constraint2).extract(time_constraint)
        u_mean = u.extract(iris.Constraint(pressure=lambda cell: 190<cell<860))
        v_mean = v.extract(iris.Constraint(pressure=lambda cell: 190<cell<860))   
        print u_mean.data
        u_mean = average_plevs(u_mean)
        v_mean = average_plevs(v_mean)
        print u_mean.data
        w = VectorWind(u_mean,v_mean)
        psi = w.streamfunction()
        psi = psi.extract(b_constraint)
        psi_list.append(psi)
    return psi_list

   

def myround(x, base=5):
    return int(base * round(float(x)/base))
    
def load_era_data(cube,dd,hr,plev,yr=2014,mth=12):
    TT = myround(hr,base=6)
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    minlat = np.amin(lats); maxlat = np.amax(lats); minlon = np.amin(lons); maxlon = np.amax(lons)
    df = '/nfs/a137/earceb/ERA-interim/{0:04d}/ap/ggap{0:04d}{1:02d}{2:02d}{3:02d}00.nc'.format(yr, mth, dd,TT)
    p_constraint = iris.Constraint(p=plev)
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    u = iris.load_cube(df,'eastward_wind')[0][:][:][:]
    v = iris.load_cube(df,'northward_wind')[0][:][:][:]
    u_mean = u.extract(iris.Constraint(p = lambda cell: 190<cell<860))
    v_mean = v.extract(iris.Constraint(p = lambda cell: 190<cell<860))    
    era_u = u_mean.collapsed('p',iris.analysis.MEAN)
    era_v = v_mean.collapsed('p',iris.analysis.MEAN)
    w = VectorWind(era_u,era_v)
    era_psi = w.streamfunction()
    era_psi = era_psi.extract(b_constraint)
    return era_psi
    
def plot_data(psi_list,avg_psi,outfile,dt,md,TT,pos_id,plev,dd,hr,dtick=10):
    #brewer_cmap = mpl_cm.get_cmap('brewer_OrRd_09')
    img_df = '/nfs/a37/scjea/storm_images/cimms/data/NWPacific/201412{0:02d}/250to850mbDeepLayerMeanLarge/201412{0:02d}.{1:02d}.NWPacific.250to850mbDeepLayerMeanLarge.png'.format(dd,hr)
    cimms_img = mpimg.imread(img_df)
    ensemble_members=np.arange(12)
    psi_levs = np.arange(-1e8,1e8,6e6)
    for em in ensemble_members:
        dataname = '/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver_glm_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(md,TT,em)
        track_data = np.load(dataname)
        if em == 0:
            lats = track_data['lats']
            lons = track_data['lons']
        else:
            lats = np.vstack([lats,track_data['lats']])
            lons = np.vstack([lons,track_data['lons']])
    plt.figure(figsize=(20,10), dpi=100)
    plt.gcf().subplots_adjust(hspace=0.05, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925)
    for em in ensemble_members:
        psi = psi_list[em]
        x = psi.coord('longitude').points
        y = psi.coord('latitude').points
        X, Y = np.meshgrid(x,y)
        plt.subplot(4,4,em+1)
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        sfunc = iplt.contour(psi,levels=psi_levs,colors='darkcyan',linewidths=2)
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
                ax.plot(lons[i,:],lats[i,:],color='lightgrey')
            else:
                mem=i
        ax.plot(lons[mem,:],lats[mem,:],color='b')
        ax.plot(lons[mem,pos_id],lats[mem,pos_id],'r*')
        map_formatter(ax,bottom_label=bottom_label, left_label=left_label, labelsize=12,tick_base_x=dtick,tick_base_y=dtick)
        ax.annotate('{0:02d}'.format(em), xy=(0.99,0.01), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
    
    plt.subplot(4,4,13)
    x = avg_psi.coord('longitude').points
    y = avg_psi.coord('latitude').points
    X, Y = np.meshgrid(x,y)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    sfunc = iplt.contour(avg_psi,levels=psi_levs,colors='darkcyan',linewidths=2)
    ax = plt.gca()
    map_formatter(ax,labelsize=12,tick_base_x=dtick,tick_base_y=dtick)
    ax.annotate('avg', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
        
    ### ERA data
    plt.subplot(4,4,16)
    plt.imshow(cimms_img)
    plt.axis('off')
    
    
    
    
    #x = era_psi.coord('longitude').points
    #y = era_psi.coord('latitude').points
    #X, Y = np.meshgrid(x,y)
    #matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    #sfunc = iplt.contour(era_psi,levels=psi_levs,colors='darkcyan',linewidths=2)
    #ax = plt.gca()
    #map_formatter(ax,labelsize=12,tick_base_x=dtick,tick_base_y=dtick,left_label=False)
    #ax.annotate('ERA', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
    
    plt.suptitle('Streamlines @ 850-200 hPa'.format(plev))    
    #time_unit, validity_time, data_time, lead_time = get_time(avg_psi)
    #ax.annotate('Simulation start time: {:%Y/%m/%d T%H:%MZ}, Lead time: +{:02d}Z'.format(data_time, int(lead_time)),xy=(0.35,0.2),horizontalalignment='left', verticalalignment='top',xycoords='figure fraction')
    ax.annotate('T+{0:03d}Z'.format(dt),xy=(0.35,0.2),horizontalalignment='left', verticalalignment='top',xycoords='figure fraction')
    
    
    
    #fig2 = plt.figure(2)
    #ax2 = fig2.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    #x = avg_psi.coord('longitude').points
    #y = avg_psi.coord('latitude').points
    #X, Y = np.meshgrid(x,y)
    #matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    #sfunc = iplt.contour(avg_psi,levels=psi_levs,colors='darkcyan',linewidths=2)
    #ax = plt.gca()
    #map_formatter(ax,labelsize=12,tick_base_x=10,tick_base_y=10)
    #ax.annotate('avg', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
    
    #psi_levs = np.arange(-1e8,1e8,2e6)
    #fig2 = plt.figure(2)
    #x = avg_psi.coord('longitude').points
    #y = avg_psi.coord('latitude').points
    #X, Y = np.meshgrid(x,y)
    #matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    #sfunc = iplt.contour(avg_psi,levels=psi_levs,colors='darkcyan',linewidths=2)
    #ax = plt.gca()
    #map_formatter(ax,labelsize=12,tick_base_x=10,tick_base_y=10)
    #ax.annotate('avg', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')

    
    
    plt.show()
    plt.savefig(outfile)
    plt.close()

def plot_avg_slines(avg_psi,outfile,md,TT,dd,hr,pos_id,dtick=10):
    ensemble_members = np.arange(12)
    psi_levs = np.arange(-1e8,1e8,2e6)
    for em in ensemble_members:
        dataname = '/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver_glm_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(md,TT,em)
        track_data = np.load(dataname)
        if em == 0:
            lats = track_data['lats']
            lons = track_data['lons']
        else:
            lats = np.vstack([lats,track_data['lats']])
            lons = np.vstack([lons,track_data['lons']])
    avg_lat = []; avg_lon = []
    for i in np.arange(len(lats[0,:])):
        avg_lat.append(np.mean(lats[:,i]))
        avg_lon.append(np.mean(lons[:,i]))
    print avg_lon
    plt.figure(figsize=(20,10), dpi=100)
    plt.gcf().subplots_adjust(hspace=0.05, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    sfunc = iplt.contour(avg_psi,levels=psi_levs,colors='darkcyan',linewidths=2)
    plt.plot(avg_lon,avg_lat,'b-',linewidth=2.5)
    plt.plot(avg_lon[pos_id],avg_lat[pos_id],'r*',MarkerSize=5)
    minlat = -10;maxlat = 30; minlon = 105; maxlon = 180;
    ax = plt.gca()
    map_formatter(ax,labelsize=35,tick_base_x=dtick,tick_base_y=dtick)
    ax.set_xlim([minlon,maxlon]); ax.set_ylim([minlat,maxlat])
    #plt.show()
    plt.savefig(outfile)
    plt.close()
    
def main():
    monthday = [1104]; times = [12]
    vt_day = [4,5,6,7,8,9]; vt_hr = [0, 6, 12, 18]
    plevs = [200]
    #### Load cubes:
    for md in monthday:
        for TT in times:
            dt = 0
            pos_id = 0
            for dd in vt_day:
                for hr in vt_hr:
                    if dd==(md-1100) and hr < TT:
                        continue
                    for plev in plevs:
                        outfile = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/plots/extended_abstract/hai_avg_streamlines_850-200hPa_{0:02d}_{1:02d}Z.png'.format(dd,hr)
                        fpath = '/nfs/a37/scjea/model_runs/Haiyan/u-ap121/data/glm/{0:04d}_{1:02d}Z/wind_plev/uvplev_glm_{0:04d}_{1:02d}Z_em'.format(md,TT)
                        print('dd = {0:02d}, hr = {1:02d}, p = {2:03d}hPa'.format(dd,hr,plev))
                        psi_list = load_ens_members(fpath,dd,hr,plev,yr=2013,mth=11)
                        avg_psi = avg_cubes(psi_list)
                        #era_psi = load_era_data(avg_psi,dd,hr,plev)
                        plot_avg_slines(avg_psi,outfile,md,TT,dd,hr,pos_id)
#                        plot_data(psi_list,avg_psi,outfile,dt,md,TT,pos_id,plev,dd,hr,dtick=15)
                    dt = dt + 6
                    pos_id = pos_id + 1


        
    
if __name__=='__main__':
    main()
