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
from vorticity import calc_div_spherical
import matplotlib.image as mpimg
from windspharm.iris import VectorWind
import numpy.ma as ma
from vorticity import calc_div_spherical

def load_ens_members(md,TT,dd,hr,mth=12,yr=2014):
    wind_plev = 200
    minlat = -2;maxlat = 55; minlon = 90; maxlon = 180;
    ens_members = np.arange(12)
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    fpath = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/glm/{0:04d}_{1:02d}Z/'.format(md,TT)
    fpath2 = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/glm/{0:04d}_{1:02d}Z/'.format(md,TT)
    olr_ext = 'gen_diags/olr_glm_{0:04d}_{1:02d}Z_'.format(md,TT)
    wind_ext = 'wind_plev/uvplev_glm_{0:04d}_{1:02d}Z_'.format(md,TT)
    #wind_ext = 'tctracker/tcver_glm_{0:04d}_{1:02d}Z_'.format(md,TT)
    olr_stash = 'm01s02i205'
    u_stash = 'm01s15i201'
    v_stash = 'm01s15i202'
    #u_stash = 'm01s03i225'
    #v_stash = 'm01s03i226'
        
    u_list = []; v_list = []; olr_list = []
    u_constraint = iris.AttributeConstraint(STASH=u_stash)
    v_constraint = iris.AttributeConstraint(STASH=v_stash)
    olr_constraint = iris.AttributeConstraint(STASH=olr_stash)
    windp_constraint = iris.Constraint(pressure=wind_plev)
    time_constraint1 = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr))
    time_constraint2 = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr+1))
    
    for em in ens_members:
        olr_df = fpath + olr_ext + 'em{0:02d}.pp'.format(em)
        wind_df = fpath2 + wind_ext + 'em{0:02d}.pp'.format(em)
        with iris.FUTURE.context(cell_datetime_objects=True):
            u = iris.load_cube(wind_df,u_constraint).extract(time_constraint1&windp_constraint)
            v = iris.load_cube(wind_df,v_constraint).extract(time_constraint1&windp_constraint)
            olr = iris.load_cube(olr_df,olr_constraint).extract(time_constraint1&b_constraint)
        u_list.append(u); v_list.append(v); olr_list.append(olr)
        #print olr.data
    return u_list, v_list, olr_list
        



def round_down(num, divisor):
    return num - (num%divisor)

def round_up(num,divisor):
    return num + divisor - (num%divisor)

def plot_data(psi_list,olr_list,divu_list,divv_list,div_list,avg_psi,avg_olr,avg_divu,avg_divv,avg_div,dd,hr,outfile,pos_id):
    step = 10; scale = 200
    print [np.amin(div_list[0].data), np.amax(div_list[0].data)]
    img_df = '/nfs/a37/scjea/storm_images/cimms/data/NWPacific/201412{0:02d}/UpperDivergenceLarge/201412{0:02d}.{1:02d}.NWPacific.UpperDivergenceLarge.png'.format(dd,hr)
    mean_olr = yr_mean(avg_olr,dd)
    cimms_img = mpimg.imread(img_df)
    
    coord_cube = divu_list[0]
    x = coord_cube.coord('longitude').points
    y = coord_cube.coord('latitude').points
    X, Y = np.meshgrid(x,y)    
    
    
    coord_cube = olr_list[0]
    olrx = coord_cube.coord('longitude').points
    olry = coord_cube.coord('latitude').points
    olrX, olrY = np.meshgrid(olrx,olry)   
    
    
    div_levs = np.arange(-100e-6,100e-6,25e-6)
    olr_cmap = mpl_cm.get_cmap('Greys')
    min_olr = 70; max_olr = 340
    
    dl = 10
    olr_levs = np.arange(min_olr,max_olr+10,10)
    maxcmap = olr_cmap(1e7); mincmap = olr_cmap(0)
    
    #brewer_cmap = mpl_cm.get_cmap('brewer_RdBu_11')
    #my_cmap = brewer_cmap([10,9,8,7,6,5,4,3,2,1,0])
    #brewer_cmap = colors.ListedColormap(my_cmap,"my_colormap")
    #maxcmap = brewer_cmap(1e7); mincmap = brewer_cmap(0)
    #olr_cmap, norm, olr_levs = normalise_cmap(min_olr,max_olr,0,dl,cmap=brewer_cmap)
    
    ens_members = np.arange(12)
    for em in ens_members:
        dataname = '/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver_glm_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(1203,12,em)
        track_data = np.load(dataname)
        if em == 0:
            lats = track_data['lats']
            lons = track_data['lons']
        else:
            lats = np.vstack([lats,track_data['lats']])
            lons = np.vstack([lons,track_data['lons']])
    
    fig=plt.figure(figsize=(10,10), dpi=100)
    plt.gcf().subplots_adjust(hspace=0.05, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925)
    olr_cmap.set_under(mincmap); olr_cmap.set_over(maxcmap)
        
    for em in ens_members:
        olr = olr_list[em]
        divu = divu_list[em]; divv = divv_list[em]
        divspeed = (divu**2 + divv**2)**0.5
        div = div_list[em]
        #print divspeed.data
        #indices = np.where(divspeed.data<2.5)
        #exit()
        #divu.data[:][:] = ma.array(divu.data[:][:],mask=False)
        #divu.data[indices] = ma.array(divu.data[indices],mask=True)
        #divv.data[:][:] = ma.array(divv.data[:][:],mask=False)
        #divv.data[indices] = ma.array(divv.data[indices],mask=True)
        #divu.data[:][:] = np.ma.masked_where(divspeed.data<25,divu.data)
        #divv.data[:][:] = np.ma.masked_where(divspeed.data<25,divv.data)
        #print divu
        #exit()
        #print olr_list[em].data[50][50]
        #print mean_olr.data[50][50]
        #olr.data[:][:] = mean_olr.data[:][:] - olr_list[em].data[:][:]
        #print olr.data[50][50]
        #print cube; print div; print divu; print divv
        if em < 2:
            em_pos = em + 1
        elif em < 4:
            em_pos = em + 3
        else:
            em_pos = em + 5
        plt.subplot(4,4,em_pos,projection=ccrs.PlateCarree())
        cube_plt = iplt.contourf(olr,cmap=olr_cmap,levels=olr_levs,extend='both')
        #div_plt = iplt.contour(div,levels=div_levs,colors='yellow',linewidths=1)
        div_arrow = plt.quiver(x[::step],y[::step],divu.data[::step,::step],divv.data[::step,::step],pivot='middle',transform=ccrs.PlateCarree(),color='y',scale=scale,linewidths=1)
        ax = plt.gca()
        if (em_pos-1) % 4 != 0:
            left_label=False
        else:
            left_label=True
        if em_pos < 13:
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
        map_formatter(ax,bottom_label=bottom_label, left_label=left_label, labelsize=12,tick_base_x=15,tick_base_y=10)
        ax.annotate('{0:02d}'.format(em), xy=(0.99,0.01), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white',fontsize=20)
    
    
    
#    plt.subplot(4,4,13,projection=ccrs.PlateCarree())
 #   olr_plt = iplt.contourf(avg_olr,levels=olr_levs,cmap=olr_cmap)
    #div_plt = iplt.contour(avg_div,levels=div_levs,colors='yellow',linewidths=1)
    #div_arrow = plt.quiver(x[::step],y[::step],avg_divu.data[::step,::step],avg_divv.data[::step,::step],pivot='middle',transform=ccrs.PlateCarree(),color='r',scale=scale,linewidths=0.25)
    
    
    plt.subplot(2,2,2)
    plt.imshow(cimms_img)
    plt.axis('off')
        
    #sfunc = iplt.contour(avg_psi,levels=psi_levs,colors='darkcyan',linewidths=2)
    #ax = plt.gca()
    #map_formatter(ax,labelsize=10,tick_base_x=10,tick_base_y=10)
    #ax.annotate('avg', xy=(1,0), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white')
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.05])
    colorbar = plt.colorbar(cube_plt, cbar_ax, orientation='horizontal')
    colorbar.ax.tick_params(labelsize=14)
    
    #plt.show()
    plt.savefig(outfile)
    plt.close()
    
def yr_mean(cube,dd):
    df = '/nfs/a37/scjea/olr_climate/olr.day.ltm.nc'
    day_of_year = (6*31) + 28 + (4*30) + dd
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    minlat = np.amin(lats); maxlat = np.amax(lats); minlon = np.amin(lons); maxlon = np.amax(lons)
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    mean_olr = iris.load_cube(df)[day_of_year].extract(b_constraint)
    target_lats = cube.coord('latitude').points
    target_lons = cube.coord('longitude').points
    sample_points = [('latitude',target_lats),('longitude', target_lons)]
    #print mean_olr
    #print cube
    mean_olr = mean_olr.interpolate(sample_points,iris.analysis.Linear())
    #print mean_olr.data
    return mean_olr
    
    

def main():
    monthday = [1203]; times = [12]
    vt_day = [3,4,5,6,7]; vt_hr = [0, 6, 12, 18]
    minlat = -2;maxlat = 55; minlon = 90; maxlon = 180;
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    
    #### Load cubes:
    for md in monthday:
        for TT in times: 
            pos_id = 1
            for dd in vt_day:
                for hr in vt_hr:
                    psi_list = []; divu_list = []; divv_list = []; div_list = []
                    if dd==(md-1200) and hr < TT:
                        continue
                    outfile = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/plots/extended_abstract/olrdivcompglm_{0:02d}_{1:02d}Z.png'.format(dd,hr)
                    u_list, v_list, olr_list = load_ens_members(md,TT,dd,hr)
                    avg_u = avg_cubes(u_list); avg_v = avg_cubes(v_list); avg_olr = avg_cubes(olr_list)
                    
                    for u,v in zip(u_list,v_list):
                        w = VectorWind(u, v)
                        divu, divv = w.irrotationalcomponent()
                        divu = divu.extract(b_constraint); divv = divv.extract(b_constraint)
                        psi = w.streamfunction()
                        psi = psi.extract(b_constraint)
                        psi_list.append(psi); divu_list.append(divu); divv_list.append(divv)
                        ured = u.extract(b_constraint); vred = v.extract(b_constraint)
                        div = calc_div_spherical(ured,vred)
                        div_list.append(div)
                    print div_list[5].data
                    avg_psi = avg_cubes(psi_list)
                    avg_divu = avg_cubes(divu_list); avg_divv = avg_cubes(divv_list)
                    avg_div = avg_cubes(div_list)
                    
                    plot_data(psi_list,olr_list,divu_list,divv_list,div_list,avg_psi,avg_olr,avg_divu,avg_divv,avg_div,dd,hr,outfile,pos_id)
                    pos_id = pos_id + 1

if __name__=='__main__':
    main()
