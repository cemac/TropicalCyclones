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
from haversine_formula import haversine_d


def load_ens_members(fpath, NN,plev, yr=2014, mth=12):
    # fpath should be everything but ens_member.pp
    ensemble_members=np.arange(12)
    #stash_code1 = 'm01s15i201'
    #stash_code2 = 'm01s15i202'
    stash_code1 = 'm01s15i300'
    stash_code2 ='m01s15i301'
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
    #time_constraint = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr))
        
    for em in ensemble_members:
        df = fpath + '{0:02d}_{1:03d}.pp'.format(em,NN)
        if plev == 0:
            with iris.FUTURE.context(cell_datetime_objects=True):
                u = iris.load_cube(df,data_constraint1).extract(time_constraint)
                v = iris.load_cube(df,data_constraint2).extract(time_constraint)
            data = iris.analysis.calculus.curl(u,v)[2]
        else:
            with iris.FUTURE.context(cell_datetime_objects=True):
                u = iris.load_cube(df,data_constraint1).extract(p_constraint&b_constraint)
                v = iris.load_cube(df,data_constraint2).extract(p_constraint&b_constraint)
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
    
def plot_data(strain,avg_strain,outfile,pos_id,plev,dtick=5):
    #brewer_cmap = mpl_cm.get_cmap('brewer_OrRd_09')
    brewer_cmap = mpl_cm.get_cmap('brewer_RdBu_11')
    md=1203; TT=12
    maxcmap = brewer_cmap(1e7); mincmap = brewer_cmap(0)
    
    ensemble_members=np.arange(12)
    #sr_min = 0; sr_max = 0
    #for em in ensemble_members:
    #    srate = strain[em]-avg_strain
    #    if np.amin(srate.data)<sr_min:
    #        sr_min = np.amin(srate.data)
    #    if np.amax(srate.data)>sr_max:
    #        sr_max=np.amax(srate.data)
    #dl = (sr_max - sr_min) / 20.
    
    sr_min = -24e-5
    sr_max = 27e-5
    dl=3e-5
    
    cmap, norm, vortlevels = normalise_cmap(sr_min,sr_max,0,dl,cmap=brewer_cmap)
    
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
    for em in ensemble_members:
        srate = strain[em] - avg_strain
        srate = strain[em]
        
        plt.subplot(3,4,em+1,projection=ccrs.PlateCarree())
        ws_plt = iplt.contourf(srate,cmap=cmap,norm=norm,levels=vortlevels,extend='both')
        ws_plt.cmap.set_under(mincmap); ws_plt.cmap.set_over(maxcmap)
        ax = plt.gca()
        if em % 4 != 0:
            left_label=False
        else:
            left_label=True
        if em < 8 :
            bottom_label=False
        else:
            bottom_label=True
        map_formatter(ax,bottom_label=bottom_label, left_label=left_label, labelsize=12,tick_base_x=10,tick_base_y=10)
        ax.annotate('{0:02d}'.format(em), xy=(0.99,0.01), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
        for i in ensemble_members:
            if i!=em:
                ax.plot(lons[i,:],lats[i,:],color='lightgrey')
            else:
                mem=i
        ax.plot(lons[mem,:],lats[mem,:],color='b')
        ax.plot(lons[mem,pos_id],lats[mem,pos_id],'r*')
    
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.05])
    colorbar = plt.colorbar(ws_plt, cbar_ax, orientation='horizontal')
    fig.suptitle('Shearing, {0:03d}hPa'.format(plev))

    
    #plt.show()
    plt.savefig(outfile)
    plt.close()

def calc_grad(data, dx):
    # Basic central differencing
    # Data should be 1-d
    grad = np.zeros(data.shape)
    grad[0] = (1/dx) * (data[1] - data[0])
    grad[-1] = (1/dx) * (data[-1] - data[-2])
    grad[1:-1] = (1/(2*dx)) * (data[2:] - data[:-2])
    return grad


def vert_derivative(cube,geop):
    # input cube and geopotential data on pressure levels - computes a very inaccurate vertical gradient.
    # Data should be of the form (pressure,latitude,longitude)
    mingeop = np.amin(geop.data)
    maxgeop = np.amax(geop.data) 
    geop_levs = np.arange(mingeop,maxgeop,100)
    
    
    p_levs = cube.coord('pressure').points
    dplev = p_levs[0:-1] - p_levs[1:]
    for plev in np.arange(len(p_levs)):
        geo_ht0 = geop.data[plev][:][:]
        geo_ht_high = geop.data[plev+1][:][:]
        geo_ht_low = geop.data[plev-1][:][:]
        diff_low = geo_ht0 - geo_ht_low
        diff_high = geo_ht_high - geo_ht0
        

def calc_stretching(u,v):
    # (1/(a*cos(lat))) * (du/dlon - d(vcos(lat))/dlat)
    
    lats = u.coord('latitude').points
    lons = u.coord('longitude').points
    LONS, LATS = np.meshgrid(lons,lats)
    
    a = 6371000. # Radius of Earth
    dlatdeg = abs(lats[1]-lats[0])
    dlat = np.pi * dlatdeg / 180. # Grid spacing (assumed constant)
    
    dlondeg = abs(lons[1]-lons[0])
    dlon = np.pi * dlondeg / 180.
    # Calculate du/dlon, i.e. for each latitude
    for i in np.arange(len(u.data)):
        dudlon = calc_grad(u.data[i,:],dlon)
        if i == 0:
            u_lon = dudlon;
        else:
            u_lon = np.vstack((u_lon,dudlon));
          
    lats_rad = LATS * np.pi / 180.
    vcoslat = v.data * np.cos(lats_rad)
    
    # Calculate d(u sin(lat))/dlat, i.e. for each lon
    for i in np.arange(len(u.data[0])):
        dvdlat = calc_grad(vcoslat[:,i],dlat)
        if i == 0:
            v_lat = dvdlat
        else:
            v_lat = np.vstack((v_lat,dvdlat))
    v_lat = np.swapaxes(v_lat,0,1)

    vrt = u[:]
    vrt.data = (1/(a*np.cos(lats_rad)))*(u_lon - v_lat)
    vrt.units = 's-1'
    vrt.rename('Atmospheric strtching')
    vrt.attributes = {'source':'Calculated using (u,v) from Met Office Unified Model', 'um_version': '10.6'}
    return vrt

def calc_shearing(u,v):
    # (1/(a*cos(lat))) * (du/dlon - d(vcos(lat))/dlat)
    
    lats = u.coord('latitude').points
    lons = u.coord('longitude').points
    LONS, LATS = np.meshgrid(lons,lats)
    
    a = 6371000. # Radius of Earth
    dlatdeg = abs(lats[1]-lats[0])
    dlat = np.pi * dlatdeg / 180. # Grid spacing (assumed constant)
    
    dlondeg = abs(lons[1]-lons[0])
    dlon = np.pi * dlondeg / 180.
    # Calculate du/dlon, i.e. for each latitude
    for i in np.arange(len(u.data)):
        dvdlon = calc_grad(v.data[i,:],dlon)
        if i == 0:
            v_lon = dvdlon;
        else:
            v_lon = np.vstack((v_lon,dvdlon));
          
    lats_rad = LATS * np.pi / 180.
    ucoslat = u.data * np.cos(lats_rad)
    
    # Calculate d(u sin(lat))/dlat, i.e. for each lon
    for i in np.arange(len(u.data[0])):
        dudlat = calc_grad(ucoslat[:,i],dlat)
        if i == 0:
            u_lat = dudlat
        else:
            u_lat = np.vstack((u_lat,dudlat))
    u_lat = np.swapaxes(u_lat,0,1)

    vrt = u[:]
    vrt.data = (1/(a*np.cos(lats_rad)))*(v_lon + u_lat)
    vrt.units = 's-1'
    vrt.rename('Atmospheric shearing')
    vrt.attributes = {'source':'Calculated using (u,v) from Met Office Unified Model', 'um_version': '10.6'}
    return vrt



def calc_strainrate(u,v):
    # Initially assuming cartesian grid
    
    lats = u.coord('latitude').points
    lons = u.coord('longitude').points
    
    LONS, LATS = np.meshgrid(lons,lats)
    
    dx = haversine_d(lats[0],lats[0],lons[0],lons[1]); 
    dy = haversine_d(lats[0],lats[1],lons[0],lons[0])
    
    # Calc derivatives
    for i in np.arange(len(u.data)):
        dvdx = calc_grad(v.data[i,:],dx)
        if i == 0:
            v_x = dvdx;
        else:
            v_x = np.vstack((v_x,dvdx));
    for i in np.arange(len(u.data)):
        dudx = calc_grad(u.data[i,:],dx)
        if i == 0:
            u_x = dudx;
        else:
            u_x = np.vstack((u_x,dudx));
            
            
    for i in np.arange(len(u.data[0])):
        dvdy = calc_grad(v.data[:,i],dy)
        if i == 0:
            v_y = dvdy;
        else:
            v_y = np.vstack((v_y,dvdy));
    for i in np.arange(len(u.data[0])):
        dudy = calc_grad(u.data[:,i],dy)
        if i == 0:
            u_y = dudy;
        else:
            u_y = np.vstack((u_y,dudy));
            
    u_y = np.swapaxes(u_y,0,1)
    v_y = np.swapaxes(v_y,0,1)
    
    strainrate = u
    strainrate.data = (2*u_x**2 + (u_y+v_x)**2 + 2*v_y**2)**0.5
    strainrate.units = 's-1'
    return strainrate
    
def main():
    monthday = [1203]; times = [12]
    vt_day = [3,4,5,6,7]; vt_hr = [0, 6, 12, 18]
    plevs = [850, 500, 200]; 
    #### Load cubes:
    for md in monthday:
        for TT in times:
            pos_id = 0
            NN=0
            for dd in vt_day:
                for hr in vt_hr:
                    NN=NN+6
                    if dd==(md-1200) and hr < TT+6:
                        continue
                    for plev in plevs:
                        strain = []
                        outfile = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/plots/glm/{0:04d}_{1:02d}Z/ensemble_plots/tcremshearing_{2:02d}_{3:02d}Z_{4:03d}hPa.png'.format(md,TT,dd,hr,plev)
                        fpath = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/glm/{0:04d}_{1:02d}Z/tcremoval/veldecomp_500km_1203_12_em'.format(md,TT)
                        u_list, v_list, min_val, max_val = load_ens_members(fpath,NN,plev)
                        for u,v in zip(u_list,v_list):
                            strain.append(calc_shearing(u,v))
                        avg_strain = avg_cubes(strain)
                        plot_data(strain,avg_strain,outfile,pos_id,plev)
                    pos_id = pos_id + 1
                    


        
    
if __name__=='__main__':
    main()
