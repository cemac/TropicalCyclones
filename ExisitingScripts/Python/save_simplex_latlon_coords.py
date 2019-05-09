import iris
import numpy as np
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from custom_cmap import *
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import datetime
from matplotlib.mlab import griddata
import matplotlib.ticker as ticker
#import Get_radial_tangent_components as rt
import colorsys
import glob
import cPickle as pickle
#######################################################
#MODEL stash data
U_CUBE_STASH='m01s15i201'  
U_FILESTREAM='b'
V_CUBE_STASH='m01s15i202'  
V_FILESTRAM='b'
SLP_CUBE_STASH='m01s16i222'  
SLP_FILESTREAM='a'
SURFACEU_CUBE_STASH='m01s03i225'
SURFACEU_FILESTREAM='a'
SURFACEV_CUBE_STASH='m01s03i226'
SURFACEV_FILESTREAM='a'
RADIATION_TRACE_STASH=['m01s00i577','Radiation (PVU)']
MICROPHYSICS_TRACE_STASH=['m01s00i580','Microphysics(PVU)']
GRAVITYWAVEDRAG_TRACE_STASH=['m01s00i581','Gravity Wave Drag (PVU) ']
SLOWPHYSICS_TRACE_STASH=['m01s00i582', 'Slow Physics (PVU)']
BOUNDARYLAYER_TRACE_STASH=['m01s00i584', 'Boundary Layer (PVU)']
CLOUD_TRACE_STASH=['m01s00i586','Cloud scheme (PVU)']
TOTALPV_TRACE_STASH=['m01s00i589', 'Total PV (PVU)']
ADVECTION_TRACE_STASH='m01s00i590'
DYNAMICALSOLVER_TRACE_STASH='m01s00i591'
MASSCONSERVATION_TRACE_STASH='m01s00i592'
INITIALPVADVECTED_TRACE_STASH=['m01s00i593', 'Initial PV advected (PVU)']
THETABOUNDARYLAYERMIXING_TRACE_STASH='m01s00i600'
THETAINTIALADVECTED_TRACE_STASH='m01s00i602'
THETABOUNDARYLAYERLATENTHEATING_TRACE_STASH='m01s00i603'
THETAMICROPHYSICS_TRACE_STASH='m01s00i605'
THETARADIATION_TRACE_STASH='m01s00i605'
MODELPV_STASH=['m01s15i218','Model PV (PVU)']
###########################################################
#Model domain data
TIMETOLOAD = 0.0   #hours after initilization use a float not an integar
INITIAL_LATITUDE=18.4
INITIAL_LONGITUDE=314.6
#INITIAL_LATITUDE=17.0
#INITIAL_LONGITUDE=133.5
DOMAIN_MIN_LON=0
DOMAIN_MAX_LON=360
DOMAIN_MIN_LAT=-90
DOMAIN_MAX_LAT=90
DOMAIN_LAT_POINTS=960
DOMAIN_LON_POINTS=1280

DOMAIN_MIN_LON=290
DOMAIN_MAX_LON=325.96
DOMAIN_MIN_LAT=8.3
DOMAIN_MAX_LAT=28.26
DOMAIN_LAT_POINTS=500
DOMAIN_LON_POINTS=900


DOMAIN_LAT_DEGREES=DOMAIN_MAX_LAT-DOMAIN_MIN_LAT
DOMAIN_LON_DEGREES=DOMAIN_MAX_LON-DOMAIN_MIN_LON
INITIAL_LATITUDE_POINTS=int(np.round((((INITIAL_LATITUDE-DOMAIN_MIN_LAT)/(DOMAIN_LAT_DEGREES+0.0))*DOMAIN_LAT_POINTS)-1))
INITIAL_LONGITUDE_POINTS=int(np.round((((INITIAL_LONGITUDE-DOMAIN_MIN_LON)/(DOMAIN_LON_DEGREES+0.0))*DOMAIN_LON_POINTS)-1))


timeconstraint=iris.Constraint(forecast_period=TIMETOLOAD)

minlonpoints=int(INITIAL_LONGITUDE_POINTS-1*(DOMAIN_LON_POINTS/(DOMAIN_LON_DEGREES+0.0)))
maxlonpoints=int(INITIAL_LONGITUDE_POINTS+1*(DOMAIN_LON_POINTS/(DOMAIN_LON_DEGREES+0.0)))
minlatpoints=int(INITIAL_LATITUDE_POINTS-1*(DOMAIN_LAT_POINTS/(DOMAIN_LAT_DEGREES+0.0)))
maxlatpoints=int(INITIAL_LATITUDE_POINTS+1*(DOMAIN_LAT_POINTS/(DOMAIN_LAT_DEGREES+0.0)))

print(minlonpoints)
print(maxlonpoints)
print(minlatpoints)
print(maxlatpoints)

longitude_constraint1=iris.Constraint(longitude=0)
#longitude_constraint1=iris.Constraint(longitude = lambda cell:cell>minlon)
longitude_constraint2=iris.Constraint(longitude = lambda cell:cell<maxlon)
latitude_constraint1=iris.Constraint(latitude = lambda cell:cell>minlat)
latitude_constraint2=iris.Constraint(latitude = lambda cell:cell<maxlat)
box_constraint=longitude_constraint1&longitude_constraint2&latitude_constraint1&latitude_constraint2



#print(1/0)



def min_vals(cube):
    # Find latitude/longitude coordinates of maximum value of data in cube.
    index = np.argmin(cube.data)
    indices= np.unravel_index(index,cube.data.shape)
    cenlat = cube.coord('latitude').points[indices[0]] # Find the coordinates for this central position
    cenlon = cube.coord('longitude').points[indices[1]]
    minval = cube.data[indices]/100.0
    return cenlat, cenlon, minval

for j in range(16,18):
   filename="wbasem{0}.pp".format(j)
   lats=[]
   lons=[]
   prs=[]
   timeconstraint=iris.Constraint(forecast_period=0)
   INITIAL_LATITUDE=18.4
   INITIAL_LONGITUDE=314.6
   INITIAL_LATITUDE_POINTS=int(np.round((((INITIAL_LATITUDE-DOMAIN_MIN_LAT)/(DOMAIN_LAT_DEGREES+0.0))*DOMAIN_LAT_POINTS)-1))
   INITIAL_LONGITUDE_POINTS=int(np.round((((INITIAL_LONGITUDE-DOMAIN_MIN_LON)/(DOMAIN_LON_DEGREES+0.0))*DOMAIN_LON_POINTS)-1))
   minlonpoints=int(INITIAL_LONGITUDE_POINTS-1*(DOMAIN_LON_POINTS/(DOMAIN_LON_DEGREES+0.0)))
   maxlonpoints=int(INITIAL_LONGITUDE_POINTS+1*(DOMAIN_LON_POINTS/(DOMAIN_LON_DEGREES+0.0)))
   minlatpoints=int(INITIAL_LATITUDE_POINTS-1*(DOMAIN_LAT_POINTS/(DOMAIN_LAT_DEGREES+0.0)))
   maxlatpoints=int(INITIAL_LATITUDE_POINTS+1*(DOMAIN_LAT_POINTS/(DOMAIN_LAT_DEGREES+0.0)))
   print(minlonpoints)
   print(maxlonpoints)
   print(minlatpoints)
   print(maxlatpoints)
   for i in range(97):
      slp_cube_constraint=iris.AttributeConstraint(STASH=SLP_CUBE_STASH)

      print(i)
      slpcube=iris.load(filename,slp_cube_constraint&timeconstraint)[0]
      #print(slpcube.coord('latitude').points)
      #print(slpcube.coord('longitude').points)
      #print(slpcube.coord('longitude').points[minlonpoints:maxlonpoints])
      #print(slpcube.coord('latitude').points[minlatpoints:maxlatpoints])
      
      slpdata=(slpcube)[minlatpoints:maxlatpoints,minlonpoints:maxlonpoints]

      c_lat, c_lon, minpress = min_vals(slpdata)
      print(c_lat)
      print(c_lon)
      print(minpress)
      #fig = plt.figure()
      #ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())

      #x = (slpdata.coord('longitude').points)-360
      #y = slpdata.coord('latitude').points
      #X,Y = np.meshgrid(x,y)
      #contour = ax.contourf(X,Y,slpdata.data)


      #gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0, color='k', linestyle=':')
      #gl.xlabels_top = False
      #gl.ylabels_right = False
      #gl.xlocator = mticker.MultipleLocator(base=2)
      #gl.ylocator = mticker.MultipleLocator(base=2)
      #gl.xlabel_style= {'size':5}
      #gl.ylabel_style= {'size':5}
      #gl.xformatter = LONGITUDE_FORMATTER
      #gl.yformatter = LATITUDE_FORMATTER
      #ax.coastlines(resolution='10m', color='k', linewidth=1)

      #plt.show()
       
      lats.append(c_lat)
      lons.append(c_lon)
      prs.append(minpress)
      INITIAL_LATITUDE_POINTS=int(np.round((((c_lat-DOMAIN_MIN_LAT)/(DOMAIN_LAT_DEGREES+0.0))*DOMAIN_LAT_POINTS)-1))
      INITIAL_LONGITUDE_POINTS=int(np.round((((c_lon-DOMAIN_MIN_LON)/(DOMAIN_LON_DEGREES+0.0))*DOMAIN_LON_POINTS)-1))
      minlonpoints=int(INITIAL_LONGITUDE_POINTS-4*(DOMAIN_LON_POINTS/(DOMAIN_LON_DEGREES+0.0)))
      maxlonpoints=int(INITIAL_LONGITUDE_POINTS+4*(DOMAIN_LON_POINTS/(DOMAIN_LON_DEGREES+0.0)))
      minlatpoints=int(INITIAL_LATITUDE_POINTS-4*(DOMAIN_LAT_POINTS/(DOMAIN_LAT_DEGREES+0.0)))
      maxlatpoints=int(INITIAL_LATITUDE_POINTS+4*(DOMAIN_LAT_POINTS/(DOMAIN_LAT_DEGREES+0.0)))
      timeconstraint=iris.Constraint(forecast_period=float(i+1))

   
   with open('latew{0:02}.pkl'.format(j),"wb") as fp:
      pickle.dump(lats,fp)

   with open('lonew{0:02}.pkl'.format(j),"wb") as fp:
      pickle.dump(lons,fp)

   with open('prsew{0:02}.pkl'.format(j),"wb") as fp:
      pickle.dump(prs,fp)
   INITIAL_LATITUDE=18.4
   INITIAL_LONGITUDE=314.6
   














