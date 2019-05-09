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

def load_ens_members(fpath, dd, hr,plev, yr=2014, mth=12):
    # fpath should be everything but ens_member.pp
    ensemble_members=np.arange(12)
    stash_code1 = 'm01s15i201'
    stash_code2 = 'm01s15i202'
    #plev = 850 #Set to 0 for single level data
    minlat = -10;maxlat = 40; minlon = 80; maxlon = 160;
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    ens_members = np.arange(12)
    min_val = 0
    max_val = 0
    u_list = []; v_list = []
    data_constraint1 = iris.AttributeConstraint(STASH=stash_code1)
    data_constraint2 = iris.AttributeConstraint(STASH=stash_code2)
    p_constraint = iris.Constraint(pressure=plev)
    time_constraint = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr))
        
    for em in ensemble_members:
        df = fpath + '{0:02d}.pp'.format(em)
        if plev == 0:
            with iris.FUTURE.context(cell_datetime_objects=True):
                u = iris.load_cube(df,data_constraint1).extract(time_constraint)
                v = iris.load_cube(df,data_constraint2).extract(time_constraint)
            data = iris.analysis.calculus.curl(u,v)[2]
        else:
            with iris.FUTURE.context(cell_datetime_objects=True):
                u = iris.load_cube(df,data_constraint1).extract(p_constraint&time_constraint&b_constraint)
                v = iris.load_cube(df,data_constraint2).extract(p_constraint&time_constraint&b_constraint)
            data = (u**2 + v**2)**0.5
        cube_max = np.amax(data.data); cube_min = np.amin(data.data)
        if cube_max>max_val:
            max_val = cube_max
        if cube_min < min_val:
            min_val = cube_min
        u_list.append(u); v_list.append(v)
    return u_list, v_list, min_val, max_val

def myround(x, base=5):
    return int(base * round(float(x)/base))
    
def load_era_data(cube,dd,hr,plev,yr=2014,mth=12):
    #plev=850
    TT = myround(hr,base=6)
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    minlat = np.amin(lats); maxlat = np.amax(lats); minlon = np.amin(lons); maxlon = np.amax(lons)
    df = '/nfs/a137/earceb/ERA-interim/{0:04d}/ap/ggap{0:04d}{1:02d}{2:02d}{3:02d}00.nc'.format(yr, mth, dd,TT)
    p_constraint = iris.Constraint(p=plev)
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    u = iris.load_cube(df,'eastward_wind').extract(p_constraint&b_constraint)[0][:][:]
    v = iris.load_cube(df,'northward_wind').extract(p_constraint&b_constraint)[0][:][:]
    return u, v
    
def plot_data(u_list,v_list,avg_u,avg_v,era_u,era_v,min_val,max_val,outfile,dtick=5):
    #brewer_cmap = mpl_cm.get_cmap('brewer_OrRd_09')
    scale=150; step=10
    cmap = mpl_cm.get_cmap('brewer_Accent_08')
    ws_levs = np.arange(0,max_val+5,5)
    
    avg_ws = (avg_u**2 + avg_v**2)**0.5
    era_ws = (era_u**2 + era_v**2)**0.5
    
    x = u_list[0].coord('longitude').points
    y = v_list[0].coord('latitude').points
    ulon = u_list[0].coord('longitude')
    transform = ulon.coord_system.as_cartopy_projection() #coordinate reference system used by data
    
    ensemble_members=np.arange(12)

    plt.figure(figsize=(20,10), dpi=100)
    plt.gcf().subplots_adjust(hspace=0.05, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925)
    for em in ensemble_members:
        u = u_list[em]; v = v_list[em]
        ws = (u**2 + v**2)**0.5
        ws.rename('windspeed')
        
        plt.subplot(4,4,em+1,projection=ccrs.PlateCarree())
        ws_plt = iplt.contourf(ws,levels=ws_levs,cmap='plasma')
        plt.quiver(x[::step],y[::step],u.data[::step,::step],v.data[::step,::step],pivot='middle',transform=transform,color='w',scale=scale)
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
        map_formatter(ax,bottom_label=bottom_label, left_label=left_label, labelsize=12,tick_base_x=10,tick_base_y=10)
        ax.annotate('{0:02d}'.format(em), xy=(0.99,0.01), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
    
    plt.subplot(4,4,13)
    ws_plt = iplt.contourf(avg_ws,levels=ws_levs,cmap='plasma')
    plt.quiver(x[::step],y[::step],avg_u.data[::step,::step],avg_v.data[::step,::step],pivot='middle',transform=transform,color='w',scale=scale)
    #ws_plt.cmap.set_under('b'); ws_plt.cmap.set_over('r')
    ax = plt.gca()
    map_formatter(ax,labelsize=12,tick_base_x=10,tick_base_y=10)
    ax.annotate('avg', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
        
    ### ERA data
    x = era_u.coord('longitude').points
    y = era_u.coord('latitude').points
    step = step 
    plt.subplot(4,4,16)
    ws_plt = iplt.contourf(era_ws,levels=ws_levs,cmap='plasma')
    plt.quiver(x[::step],y[::step],era_u.data[::step,::step],era_v.data[::step,::step],pivot='middle',transform=transform,color='w',scale=scale)
#    ws_plt.cmap.set_under('b'); ws_plt.cmap.set_over('r')
    ax = plt.gca()
    map_formatter(ax,labelsize=12,tick_base_x=10,tick_base_y=10,left_label=False)
    ax.annotate('ERA', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
    
    colorbar_axes = plt.gcf().add_axes([0.35, 0.05, 0.3, 0.05])
    colorbar = plt.colorbar(ws_plt, colorbar_axes, orientation='horizontal')
    plt.suptitle('Windspeed')    
    #time_unit, validity_time, data_time, lead_time = get_time(avg_data)
    #ax.annotate('Simulation start time: {:%Y/%m/%d T%H:%MZ}, Lead time: +{:02d}Z'.format(data_time, int(lead_time)),xy=(0.35,0.2),horizontalalignment='left', verticalalignment='top',xycoords='figure fraction')
    #iplt.show()
    plt.savefig(outfile)
    plt.close()
    
    
def main():
    monthday = [1203]; times = [12]
    vt_day = [3,4,5,6,7]; vt_hr = [0, 6, 12, 18]
    plevs = [850, 500, 200]
    #### Load cubes:
    for md in monthday:
        for TT in times:
            for dd in vt_day:
                for hr in vt_hr:
                    if dd==(md-1200) and hr < TT:
                        continue
                    for plev in plevs:
                        outfile = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/plots/glm/{0:04d}_{1:02d}Z/ensemble_plots/winds_{2:02d}_{3:02d}Z_{4:03d}hPa.png'.format(md,TT,dd,hr,plev)
                        fpath = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/glm/{0:04d}_{1:02d}Z/wind_plev/uvplev_glm_{0:04d}_{1:02d}Z_em'.format(md,TT)
                        u_list, v_list, min_val, max_val = load_ens_members(fpath,dd,hr,plev)
                        avg_u = avg_cubes(u_list); avg_v = avg_cubes(v_list)
                        era_u, era_v = load_era_data(avg_u,dd,hr,plev)
                        plot_data(u_list,v_list,avg_u,avg_v,era_u,era_v,min_val,max_val,outfile)
                    
                    


        
    
if __name__=='__main__':
    main()
