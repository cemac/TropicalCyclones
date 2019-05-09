# Calculates the vorticity or divergence using cartesian grid
import os
import matplotlib
#matplotlib.use('Agg')
import iris
import iris.analysis.calculus
import numpy as np
from cube_tools import *
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from plotting_tools import *
from poisson2 import poisson
#from poisson2 import poisson as poisson2


def calc_grad(data, dx):
    # Basic central differencing
    # Data should be 1-d
    grad = np.zeros(data.shape)
    grad[0] = (1/dx) * (data[1] - data[0])
    grad[-1] = (1/dx) * (data[-1] - data[-2])
    grad[1:-1] = (1/(2*dx)) * (data[2:] - data[:-2])
    return grad

def plot_vort(vrt,u,v,outfile,PSI):

    psi = poisson(vrt,u,v)

    x = psi.coord('longitude').points
    y = psi.coord('latitude').points
    X,Y = np.meshgrid(x,y)

    x2 = PSI.coord('longitude').points
    y2 = PSI.coord('latitude').points
    X2,Y2 = np.meshgrid(x2,y2)

    #vrtmax = 60e-5
    #vrtmin = -60e-5
    vrtmax = 30e-5
    vrtmin = -30e-5
    dl = 6e-5
    cmap, norm, vortlevels = normalise_cmap(vrtmin,vrtmax,0,dl)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.75, color='k', linestyle=':')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlocator = mticker.MultipleLocator(base=15.0)
    gl.ylocator = mticker.MultipleLocator(base=15.0)
    gl.xlabel_style= {'size':15}
    gl.ylabel_style= {'size':15}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.coastlines(resolution='10m', color='k', linewidth=1)

    psi_levs = np.arange(-1e3,1e3,6e1)

    vort = ax.contour(X,Y,psi.data,colors='k')

    vort = ax.contour(X2,Y2,PSI.data,colors='r')


#    vort.cmap.set_under('b'); vort.cmap.set_over('r')
#    cbar = plt.colorbar(vort,orientation='vertical',fraction=0.046, pad=0.09,format=OOMFormatter(-5, mathText=False))
#    cbar.set_label('Relative vorticity, (s$^{-1}$)')

    plt.show()
    #plt.savefig(outfile)
    plt.close()
    return vort

def plot_divergence(div,outfile):

    x = div.coord('longitude').points
    y = div.coord('latitude').points
    X,Y = np.meshgrid(x,y)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.75, color='k', linestyle=':')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlocator = mticker.MultipleLocator(base=15.0)
    gl.ylocator = mticker.MultipleLocator(base=15.0)
    gl.xlabel_style= {'size':15}
    gl.ylabel_style= {'size':15}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.coastlines(resolution='10m', color='k', linewidth=1)
    div_plt = ax.pcolormesh(X,Y,div.data,vmin=-25e-4,vmax=25e-4)
    cbar = plt.colorbar(div_plt,orientation='vertical',fraction=0.046, pad=0.09,format=OOMFormatter(-4, mathText=False))

    plt.show()
    plt.savefig(outfile)
    plt.close()

def calc_div_spherical(u,v):
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
    vrt.data = (1/(a*np.cos(lats_rad)))*(u_lon + v_lat)
    vrt.units = 's-1'
    vrt.rename('Atmospheric divergence')
    vrt.attributes = {'source':'Calculated using (u,v) from Met Office Unified Model', 'um_version': '10.6'}
    return vrt


def calc_vrt_spherical(u,v):
    # Calculates vorticity in spherical coordinates using the following:
    # omega = (1/(a cos(lat)))*(dv/dlon - d(u cos(lat))/dlat)
    # u, v cubes on one vertical level

    lats = u.coord('latitude').points
    lons = u.coord('longitude').points
    LONS, LATS = np.meshgrid(lons, lats)

    a = 6371000. # Radius of Earth
    ddeg = abs(lats[1]-lats[0])
    drad = np.pi * ddeg / 180. # Grid spacing (assumed constant)

    # Calculate dv/dlon, i.e. for each latitude
    for i in np.arange(len(u.data)):
        dvdlon = calc_grad(v.data[i,:],drad)
        if i == 0:
            v_lon = dvdlon;
        else:
            v_lon = np.vstack((v_lon,dvdlon));

    lats_rad = LATS * np.pi / 180.
    ucoslat = u.data * np.cos(lats_rad)

    # Calculate d(u sin(lat))/dlat, i.e. for each lon
    for i in np.arange(len(u.data[0])):
        dudlat = calc_grad(ucoslat[:,i],drad)
        if i == 0:
            u_lat = dudlat
        else:
            u_lat = np.vstack((u_lat,dudlat))
    u_lat = np.swapaxes(u_lat,0,1)

    vrt = u[:]
    vrt.data = (1/(a*np.cos(lats_rad)))*(v_lon - u_lat)
    vrt.units = 's-1'
    vrt.standard_name='atmosphere_relative_vorticity'
    vrt.long_name = 'Calculated using spherical coordinates finite difference'
    vrt.attributes = {'source':'Calculated using (u,v) from Met Office Unified Model', 'um_version': '10.6'}
    return vrt


def main():
    plev = 200
    lead_times = np.arange(4)*3

# Load all data
    filename = "/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/4p4/1203_12Z/tcremoval/veldecomp_300km_4p4_1203_12_em05_048.pp"
    #constraint = iris.AttributeConstraint(STASH='m01s15i201')
    minlat = -5;maxlat = 37; minlon = 100; maxlon = 160;
    b_constraint=box_constraint(minlat,maxlat,minlon,maxlon)
    constraint = iris.AttributeConstraint(STASH='m01s15i304')
    cubes = iris.load(filename)
    #print cubes
    #print cubes[0]
    #print cubes[1]
    #exit()
    uwind = iris.load_cube(filename,constraint).extract(b_constraint)
    #uwind = uwind[4][:][:][:]
    #constraint = iris.AttributeConstraint(STASH='m01s15i202')
    constraint = iris.AttributeConstraint(STASH='m01s15i305')
    vwind = iris.load_cube(filename,constraint).extract(b_constraint)
    #vwind = vwind[4][:][:][:]
    uwind = uwind.extract(iris.Constraint(pressure=plev))
    vwind = vwind.extract(iris.Constraint(pressure=plev))
    psi_constraint = iris.AttributeConstraint(STASH='m01s15i307')
    psi = iris.load_cube(filename,psi_constraint).extract(b_constraint&iris.Constraint(pressure=plev))
    print psi.data
    vrt = calc_vrt_spherical(uwind,vwind)
    outfile='dud'
    plot_vort(vrt,uwind,vwind,outfile,psi)



if __name__=='__main__':
	main()
