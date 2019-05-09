# Calculates the vorticity or divergence using cartesian grid
import os
import matplotlib
#matplotlib.use('Agg')
import iris
import iris.analysis.calculus
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from custom_cmap import *
import iris.plot as qplt
import cPickle as pickle
import matplotlib.image as image


def calculate_pot_temp(pressure,temperature):
    return temperature*(1000/pressure)**(0.286)
    
    
def calc_grad(data, dx):
    # Basic central differencing
    # Data should be 1-d
    grad = np.zeros(data.shape)
    grad[0] = (1/dx) * (data[1] - data[0])
    grad[-1] = (1/dx) * (data[-1] - data[-2])
    grad[1:-1] = (1/(2*dx)) * (data[2:] - data[:-2])
    return grad


def main():
    plev = 500
    timeselect=86 #hours from start of run
    levelsarray=100*(((np.arange(8))*16)+900)
    lead_times = np.arange(2)*3
    filelist=[] #create list of files
    numfilestoload=20
    ensemble=0
    day=19
    time=0
   
    
    

    fig = plt.figure(figsize=(5,4))
    levelsvv=((np.arange(25))*10)+100

    dataarray=np.zeros((750,1000))
    meandataarray=np.zeros((250,249))
    for i in range(18):
       with open('late{0:02}.pkl'.format(i),"rb") as fp:
               latt=np.asarray(pickle.load(fp))
    
       with open('lone{0:02}.pkl'.format(i),"rb") as fp:
               lonn=np.asarray(pickle.load(fp))
       
    
       #lonn=lonn-360
       filename = "olralle{0}.pp".format(i)
       #filename2 = "/nfs/a299/TCs/maria/MARIA_09{1:02}_{2:02}Z_em{0:02}_pb.pp".format(i,day,time)
    
       #print(iris.load(filename))
       #print("pa")
       #print(iris.load(filename2))
       #print("pb")
    #for i in range(numfilestoload-1):

       ir=iris.load(filename)[0][3]
       
      
       
       minlon=float(lonn[timeselect]-5)
       maxlon=float(lonn[timeselect]+5)
       minlat=float(latt[timeselect]-5)
       maxlat=float(latt[timeselect]+5)

       print(minlon)
       print(maxlon)
       print(minlat)
       print(maxlat)

       if i==0:

          irtemplate=ir.extract(iris.Constraint(longitude = lambda cell:minlon<cell<maxlon,latitude = lambda cell:minlat<cell<maxlat))
          ir=ir.extract(iris.Constraint(longitude = lambda cell:minlon<cell<maxlon,latitude = lambda cell:minlat<cell<maxlat))

       else:
          ir=ir.extract(iris.Constraint(longitude = lambda cell:minlon<cell<maxlon,latitude = lambda cell:minlat<cell<maxlat))

       print(ir)
       #print(ir)
       #print(ir)
       #vv=vv.extract(iris.Constraint(pressure=plev))
       dataarray=ir.data
       
    #/nfs/a299/TCs/maria/MARIA_09{1:02}_{2:02}Z_em{0:02}_pa.pp
    #filename = "/nfs/a37/scjea/for_will/20170918T0000Z_MARIA_km4p4_ra1t_em00_pa{0:03}.pp".format(timeselect)
       
    
       #slp=iris.load(filename)[13] #13th cube happens to be SLP
     
       x = ir.coord('longitude').points
       y = ir.coord('latitude').points
       #print(x)
       #print(y)
       X,Y = np.meshgrid(x,y)
       #arraytimeselect=(timeselect/6)
    #sstmax = 32
    #sstmin = 26.5
    #dl = 0.5
    #cmap, norm, sstlevels = normalise_cmap(sstmin,sstmax,0,dl)
    #levelssst=[26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32]
       color_map  =  plt.cm.jet

       print(i)
       ax = fig.add_subplot(5,4,i+1, projection=ccrs.PlateCarree())
    
       gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0, color='k', linestyle=':')
       gl.xlabels_top = False
       gl.ylabels_right = False
       gl.xlocator = mticker.MultipleLocator(base=2)
       gl.ylocator = mticker.MultipleLocator(base=2)
       gl.xlabel_style= {'size':5}
       gl.ylabel_style= {'size':5}
       gl.xformatter = LONGITUDE_FORMATTER
       gl.yformatter = LATITUDE_FORMATTER
       ax.coastlines(resolution='10m', color='k', linewidth=1)
       ax.set_title("P{0}".format(i))
   
    
    #ax.set_extent([-80.0, -50.0, 10.0, 35.0]) 
       ax.set_extent([lonn[timeselect]-5, lonn[timeselect]+5, latt[timeselect]-5, latt[timeselect]+5]) 
       
    #sstcontour = ax.contourf(X,Y,sst.data,levels=sstlevels,norm=norm,cmap=cmap,extend='both')
      
       #print(1/0)
       print(dataarray.shape)
       vvcontour = ax.contourf(X,Y,dataarray,levels=levelsvv,cmap=plt.cm.binary,extend='both')
       cbar = plt.colorbar(vvcontour)
    
       #contour=qplt.contour(slp[timeselect],levelsarray,linewidths=1.0,colors='black',linestyles='-')
       newir=ir.regrid(irtemplate,iris.analysis.Linear())
       print(newir)
    
       #meandataarray=meandataarray+newir.data
       #print(meandataarray)
    
    x2 = newir.coord('longitude').points
    y2 = newir.coord('latitude').points
    #print(newir)
    
    #print(x)
    #print(y)
    X2,Y2 = np.meshgrid(x2,y2)
    
    #ax=fig.add_subplot(5,4,18)
    #vvcontour = ax.contourf(X2,Y2,meandataarray/(i+1),levels=levelsvv,cmap=plt.cm.binary,extend='both')
    #cbar = plt.colorbar(vvcontour)
    #ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)
    #ax.set_title("Storm relative mean",fontsize=8)


    ax=fig.add_subplot(5,4,19)
    ax.imshow(image.imread("20170906IRMA.png"))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_title("Visible Satellite",fontsize=8)
    #gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.75, color='k', linestyle=':')
    #gl.xlabels_top = False
    #gl.ylabels_right = False
    #gl.xlocator = mticker.MultipleLocator(base=1)
    #gl.ylocator = mticker.MultipleLocator(base=1)
    #gl.xlabel_style= {'size':10}
    #gl.ylabel_style= {'size':10}
    #gl.xformatter = LONGITUDE_FORMATTER
    #gl.yformatter = LATITUDE_FORMATTER
    #ax.coastlines(resolution='10m', color='k', linewidth=1)
    #plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.tight_layout(h_pad=0.2)
    plt.suptitle("Outgoing long wave radiation ($W/m^2$) at {0}:00".format(timeselect,plev))
    plt.show()
if __name__=='__main__':
	main()
