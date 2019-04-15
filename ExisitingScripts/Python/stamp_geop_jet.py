import iris
import iris.analysis.calculus
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from plotting_tools import *
from cube_tools import *
import matplotlib.cm as mpl_cm
import iris.plot as iplt

def load_ens_members(wind_fpath, geop_fpath, dd, hr, yr=2014, mth=12):
    # fpath should be everything but ens_member.pp
    minlat = -5;maxlat = 40; minlon = 80; maxlon = 160;
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    plev = 100
    ensemble_members=np.arange(12)
    u_stash = 'm01s15i201'
    v_stash = 'm01s15i202'
    geop_stash = 'm01s16i202'
    geop_plev = 500 #Set to 0 for single level data
    wind_plev = 200
    ens_members = np.arange(12)
    min_val = 0
    max_val = 0
    u_list = []; v_list = []; geop_list = []
    u_constraint = iris.AttributeConstraint(STASH=u_stash)
    v_constraint = iris.AttributeConstraint(STASH=v_stash)
    geop_constraint = iris.AttributeConstraint(STASH=geop_stash)
    geopp_constraint = iris.Constraint(pressure=geop_plev)
    windp_constraint = iris.Constraint(pressure=wind_plev)
    time_constraint = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr))
    time_constraintA = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr))
        
    for em in ensemble_members:
        wind_df = wind_fpath + '{0:02d}.pp'.format(em)
        geop_df = geop_fpath + '{0:02d}.pp'.format(em)
        if plev == 0:
            with iris.FUTURE.context(cell_datetime_objects=True):
                u = iris.load_cube(df,data_constraint).extract(time_constraint)
                v = iris.load_cube(df,data_constraint).extract(time_constraint)
            data = iris.analysis.calculus.curl(u,v)[2]
        else:
            with iris.FUTURE.context(cell_datetime_objects=True):
                u = iris.load_cube(wind_df,u_constraint).extract(windp_constraint&time_constraint&b_constraint)
                v = iris.load_cube(wind_df,v_constraint).extract(windp_constraint&time_constraint&b_constraint)
                geop = iris.load_cube(geop_df,geop_constraint).extract(geopp_constraint&time_constraintA&b_constraint)
                print geop
        u_list.append(u); v_list.append(v); geop_list.append(geop)
    return u_list, v_list, geop_list

def myround(x, base=5):
    return int(base * round(float(x)/base))
    
def load_era_data(cube,dd,hr,yr=2014,mth=12):
    geop_plev = 500; wind_plev = 200
    TT = myround(hr,base=6)
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    minlat = np.amin(lats); maxlat = np.amax(lats); minlon = np.amin(lons); maxlon = np.amax(lons)
    df = '/nfs/a137/earceb/ERA-interim/{0:04d}/ap/ggap{0:04d}{1:02d}{2:02d}{3:02d}00.nc'.format(yr, mth, dd,TT)
    geopp_constraint = iris.Constraint(p=geop_plev)
    windp_constraint = iris.Constraint(p=wind_plev)
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    u = iris.load_cube(df,'eastward_wind').extract(windp_constraint&b_constraint)[0][:][:]
    v = iris.load_cube(df,'northward_wind').extract(windp_constraint&b_constraint)[0][:][:]
    geop = iris.load_cube(df,'geopotential').extract(geopp_constraint&b_constraint)[0][:][:]
    return u,v,geop
    
def plot_data(u_list,v_list,geop_list,era_u,era_v,era_geop,outfile,dtick=10):
    #brewer_cmap = mpl_cm.get_cmap('brewer_OrRd_09')
    print geop_list[0]
    print era_geop
    era_geop.data = era_geop.data / 9.81
    avg_u = avg_cubes(u_list)
    avg_v = avg_cubes(v_list)
    avg_geop = avg_cubes(geop_list)
    
    ensemble_members=np.arange(12) 

    plt.figure(figsize=(20,10), dpi=100)
    plt.gcf().subplots_adjust(hspace=0.05, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925)
    for em in ensemble_members:
        plt.subplot(4,4,em+1)
        geop = geop_list[em]        
        mingeop = np.floor_divide(np.amin(geop.data),40) * 40
        maxgeop = (np.floor_divide(np.amax(geop.data),40) + 2) * 40
        geoplevels = np.arange(mingeop,maxgeop,40)
        geop = iplt.contour(geop,levels=geoplevels,colors='k')
        plt.clabel(geop)
        u = u_list[em]
        v = v_list[em]
        ws = (u**2 + v**2)**0.5
        ws.data = np.ma.masked_less(ws.data,40)
        ws = iplt.contourf(ws)
        #ws.cmap.set_under('b'); ws.cmap.set_over('r')
        ax = plt.gca()
        if em % 4 != 0:
            left_label=False
        else:
            left_label=True
        if em < 9 or em == 11:
            bottom_label=False
        else:
            bottom_label=True
        map_formatter(ax,bottom_label=bottom_label, left_label=left_label, labelsize=12,tick_base_x=dtick,tick_base_y=dtick)
        ax.annotate('{0:02d}'.format(em), xy=(0.99,0.01), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
    
    plt.subplot(4,4,13)
    mingeop = np.floor_divide(np.amin(avg_geop.data),40) * 40
    maxgeop = (np.floor_divide(np.amax(avg_geop.data),40) + 2) * 40
    geoplevels = np.arange(mingeop,maxgeop,40)
    geop = iplt.contour(avg_geop,levels=geoplevels,colors='k')
    plt.clabel(geop)
    ws = (avg_u**2 + avg_v**2)**0.5
    ws.data = np.ma.masked_less(ws.data,20)
    ws = iplt.contourf(ws)
    ax = plt.gca()
    map_formatter(ax,labelsize=12,tick_base_x=10,tick_base_y=10)
    ax.annotate('avg', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
    #    
    #### ERA data
    plt.subplot(4,4,16)
    mingeop = np.floor_divide(np.amin(era_geop.data),40) * 40
    maxgeop = (np.floor_divide(np.amax(era_geop.data),40) + 2) * 40
    geoplevels = np.arange(mingeop,maxgeop,40)
    geop = iplt.contour(era_geop,levels=geoplevels,colors='k')
    plt.clabel(geop)
    ws = (era_u**2 + era_v**2)**0.5
    ws.data = np.ma.masked_less(ws.data,20)
    ws = iplt.contourf(ws)
    ax = plt.gca()
    map_formatter(ax,labelsize=12,tick_base_x=10,tick_base_y=10,left_label=False)
    ax.annotate('ERA', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
   # 
    colorbar_axes = plt.gcf().add_axes([0.35, 0.05, 0.3, 0.05])
    colorbar = plt.colorbar(ws, colorbar_axes, orientation='horizontal')
    plt.suptitle('Windspeed > 40m s-1, geop @ 500hPa')    
    #time_unit, validity_time, data_time, lead_time = get_time(avg_data)
    #ax.annotate('Simulation start time: {:%Y/%m/%d T%H:%MZ}, Lead time: +{:02d}Z'.format(data_time, int(lead_time)),xy=(0.35,0.2),horizontalalignment='left', verticalalignment='top',xycoords='figure fraction')
    #iplt.show()
    plt.savefig(outfile)
    plt.close()
    
    
def main():
    monthday = [1203]; times = [12,0]
    vt_day = [3,4,5,6,7]; vt_hr = [0, 6, 12, 18]
    #### Load cubes:
    for md in monthday:
        for TT in times:
            for dd in vt_day:
                for hr in vt_hr:
                    if dd==(md-1200) and hr < TT:
                        continue
                    outfile = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/plots/glm/{0:04d}_{1:02d}Z/ensemble_plots/geopjet_500hPa_{2:02d}_{3:02d}Z.png'.format(md,TT,dd,hr)
                    wind_fpath = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/glm/{0:04d}_{1:02d}Z/wind_plev/uvplev_glm_{0:04d}_{1:02d}Z_em'.format(md,TT)
                    geop_fpath = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/glm/{0:04d}_{1:02d}Z/geop/geop_glm_{0:04d}_{1:02d}Z_em'.format(md,TT)
                    u_list, v_list, geop_list = load_ens_members(wind_fpath,geop_fpath,dd,hr)
                    era_u, era_v, era_geop = load_era_data(u_list[0],dd,hr)
                    plot_data(u_list,v_list,geop_list,era_u,era_v,era_geop,outfile,dtick=15)
                    
                    


        
    
if __name__=='__main__':
    main()
