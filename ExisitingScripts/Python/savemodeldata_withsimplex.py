

################################
import iris
import numpy as np
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from custom_cmap import *
import cartopy.crs as ccrs

import matplotlib.ticker as mticker

###############################

#FILE data

FOLDER_LOCATION="/nfs/a299/TCs/neptark/tracers/phi_4p4_20160703T1200Z" #wherethefilesarestored

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
RADIATION_TRACE_STASH='m01s00i577'
MICROPHYSICS_TRACE_STASH='m01s00i580'
GRAVITYWAVEDRAG_TRACE_STASH='m01s00i581'
SLOWPHYSICS_TRACE_STASH='m01s00i582'
BOUNDARYLAYER_TRACE_STASH='m01s00i584'
CLOUD_TRACE_STASH=['m01s00i586','Cloud scheme']
TOTALPV_TRACE_STASH='m01s00i589'
ADVECTION_TRACE_STASH='m01s00i590'
DYNAMICALSOLVER_TRACE_STASH='m01s00i591'
MASSCONSERVATION_TRACE_STASH='m01s00i592'
INITIALPVADVECTED_TRACE_STASH='m01s00i593'
THETABOUNDARYLAYERMIXING_TRACE_STASH='m01s00i600'
THETAINTIALADVECTED_TRACE_STASH='m01s00i602'
THETABOUNDARYLAYERLATENTHEATING_TRACE_STASH='m01s00i603'
THETAMICROPHYSICS_TRACE_STASH='m01s00i605'
THETARADIATION_TRACE_STASH='m01s00i605'

#Model domain data
#TIMETOLOAD = 93.0   #hours after initilization use a float not an integar
INITIAL_LATITUDE=17.2
INITIAL_LONGITUDE=133.5

import cPickle as pickle




def cylindrical(stash_name,xcoordinate,ycoordinate,radius=200,ztickspacing=100,rtickspacing=20):
    phi_interval = np.pi / 16  #azimuthal intervals
    rad_interval=10.0
    ranges = np.arange(0,radius,rad_interval)   # from 0 to 500km with 10km intervals
    phis = np.arange(0,2*np.pi,phi_interval) #phi angles, not sure about theta
    [x0, y0] = xcoordinate,ycoordinate
    
    
    cube_constraint=iris.AttributeConstraint(STASH=stash_name[0])
    variable=cubeattime.extract(cube_constraint)[0]

    try:
      z=variable.coord('model_level_number').points
      height=variable.coord('level_height').points
    except iris.exceptions.CoordinateNotFoundError: 
       pass
    z=variable.coord('model_level_number').points
    
    height=variable.coord('level_height').points
    
    #print(1/0)
    final=np.zeros((z.size,ranges.size))
    for levels in range(z.size):
      
       functofZ_variable=variable.extract(iris.Constraint(model_level_number=z[levels]))
       #print(functofZ_variable.data)
       for r in ranges:
           for phi in phis:
               xpoi = x0 + (1/(111.32*np.cos((y0*np.pi)/180.0)))*r*np.cos(phi) #Check the 0.009  0.009 degrees is 1km
               #xpoi = x0 + 0.009*r*np.cos(phi)
               ypoi = y0 + 0.009*r*np.sin(phi)
           
               new_point = [('latitude',ypoi),('longitude',xpoi)]
               
               functofZRPhi=functofZ_variable.interpolate(new_point,iris.analysis.Linear()).data
               #print(functofZRPhi)
          
               if phi == 0:    #I'm guessing this sets it so 0 degrees is always in the direction of the storm motion. 
                  
                   functofZR=functofZRPhi
               else:
                   
                   functofZR=np.append(functofZR,functofZRPhi)
                   
        # Average everything in this section
           if r == 0:
              
               functofZ=np.mean(functofZR)
           else:
              
               functofZ=np.append(functofZ,np.mean(functofZR))
           print(r)

       final[levels]=functofZ
       #print(final[levels])
    fig = plt.figure()
    
    ax = fig.add_subplot(1,1,1)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(ztickspacing))
    ax.xaxis.set_major_locator(mticker.MultipleLocator(rtickspacing))
    plt.xticks(rotation=45)
    R,Z = np.meshgrid(ranges,height)
   
    term1=ax.contourf(R,Z,final)
    ax.set_title(stash_name[1])
    ax.set_xlabel('Radial distance from storm centre (km)')
    ax.set_ylabel('Hybrid height level (m)')
    cbar = plt.colorbar(term1)
    plt.show()

def cylindrical_azirad(variable_u,variable_v,xcoordinate,ycoordinate,radius=200,ztickspacing=1000,rtickspacing=20):
    phi_interval = np.pi / 16  #azimuthal intervals
    rad_interval=10.0
    ranges = np.arange(0,radius,rad_interval)   # from 0 to 500km with 10km intervals
    phis = np.arange(0,2*np.pi,phi_interval) #phi angles, not sure about theta
    [x0, y0] = xcoordinate,ycoordinate
  
    #cube_constraint_u=iris.AttributeConstraint(STASH=stash_u_name[0])
    #cube_constraint_v=iris.AttributeConstraint(STASH=stash_v_name[0])
    #variable_u=cubeattime.extract(cube_constraint_u)[0]
    #variable_v=cubeattime.extract(cube_constraint_v)[0]
    try:
      z=variable_u.coord('model_level_number').points
      
      
    except iris.exceptions.CoordinateNotFoundError: 
       try:
          z=variable_u.coord('level_height').points
       except iris.exceptions.CoordinateNotFoundError: 
          z=variable_u.coord('pressure').points



    
    azi_final=np.zeros((z.size,ranges.size))
    rad_final=np.zeros((z.size,ranges.size))
    for levels in range(z.size):
      
       functofZ_variable_u=(variable_u)[levels]
       functofZ_variable_v=(variable_v)[levels]
  
       for r in ranges:
           for phi in phis:
               xpoi = x0 + (1/(111.32*np.cos((y0*np.pi)/180.0)))*r*np.cos(phi) #Check the 0.009  0.009 degrees is 1km
               #xpoi = x0 + 0.009*r*np.cos(phi)
               ypoi = y0 + 0.009*r*np.sin(phi)
           
               new_point = [('latitude',ypoi),('longitude',xpoi)]
               functofZRPhi_u=functofZ_variable_u.interpolate(new_point,iris.analysis.Linear()).data
               functofZRPhi_v=functofZ_variable_v.interpolate(new_point,iris.analysis.Linear()).data
               azi_functofZRPhi = -functofZRPhi_u*np.sin(phi) + functofZRPhi_v*np.cos(phi)
               rad_functofZRPhi =  functofZRPhi_u*np.cos(phi) + functofZRPhi_v*np.sin(phi)
               
              
               #print(functofZRPhi)
          
               if phi == 0:    #I'm guessing this sets it so 0 degrees is always in the direction of the storm motion. 
                  
                   azi_functofZR=azi_functofZRPhi
                   rad_functofZR=rad_functofZRPhi
               else:
                   
                   azi_functofZR=np.append(azi_functofZR,azi_functofZRPhi)
                   rad_functofZR=np.append(rad_functofZR,rad_functofZRPhi)
                   
        # Average everything in this section
           if r == 0:
              
               azi_functofZ=np.mean(azi_functofZR)
               rad_functofZ=np.mean(rad_functofZR)
           else:
              
               azi_functofZ=np.append(azi_functofZ,np.mean(azi_functofZR))
               rad_functofZ=np.append(rad_functofZ,np.mean(rad_functofZR))
           

       azi_final[levels]=azi_functofZ
       rad_final[levels]=rad_functofZ
       #print(final[levels])
    
    return ranges,z,azi_final,rad_final
   






u_cube_constraint=iris.AttributeConstraint(STASH=U_CUBE_STASH)
v_cube_constraint=iris.AttributeConstraint(STASH=V_CUBE_STASH)

#fig,axes = plt.subplots(nrows=5,ncols=4,sharex=True,sharey=True)
#plt.setp(plt.xticks()[1],rotation=45)
e=0
t=36

#for ax in axes.flat:
#for e in range(18):
   #for t in range(96):   
print(e)
print(t)
   #if e<18: 
TIMETOLOAD=t
with open('/nfs/a319/ee16wst/model_track/late{0:02}.pkl'.format(e),"rb") as fp:
   latt=np.asarray(pickle.load(fp))
    
with open('/nfs/a319/ee16wst/model_track/lone{0:02}.pkl'.format(e),"rb") as fp:
   lonn=np.asarray(pickle.load(fp))
   
#with open('/nfs/a319/ee16wst/late_simplex{0:02}.pkl'.format(e),"rb") as fp:
   #latt=np.asarray(pickle.load(fp))
    
#with open('/nfs/a319/ee16wst/lone_simplex{0:02}.pkl'.format(e),"rb") as fp:
   #lonn=np.asarray(pickle.load(fp))
   
#print(lonn[36])
#print(latt[36])

#print("HERE")
#print lonn[36]
#print("HERE")
#print(1/0)
filename="/nfs/a319/ee16wst/model_data/basem{0}.pp".format(e)
timeconstraint=iris.Constraint(forecast_period=TIMETOLOAD)
ucube=iris.load(filename,u_cube_constraint)[0].extract(timeconstraint) 
vcube=iris.load(filename,v_cube_constraint)[0].extract(timeconstraint) 
ranges,z,azi_final,rad_final=cylindrical_azirad(ucube,vcube,360-52.7012941957,16.2955161746)
np.save("cyl_azi_e{0}_t{1}_sx.npy".format(e,t),azi_final)
np.save("cyl_rad_e{0}_t{1}_sx.npy".format(e,t),rad_final)
      
      #a=np.load("aziradtest.npy")
      #print(a)
      #print("gap")
      #print(azi_final)
      #print(1/0)
      #R,Z = np.meshgrid(ranges,z)
   
      #term1=ax.contourf(R,Z,azi_final,cmap=plt.cm.jet,levels=[-10,0,10,20,30,40,50,60,70,80,90],extend='both')
    
   
   #ax.yaxis.set_major_locator(mticker.MultipleLocator(100))
   #ax.xaxis.set_major_locator(mticker.MultipleLocator(20))
   
  
   
   #ax.set_title(stash_name[1])
   
   #ax.annotate('P{0}'.format(e), xy=(0.97,0.03), xycoords='axes fraction',horizontalalignment='right',verticalalignment='bottom',color='k',backgroundcolor='white',fontsize=20)
   #e+=1
   

#ax.invert_yaxis()   
#ax.set_xlabel('Radial distance from storm centre (km)')
#ax.set_ylabel('Pressure level (hPa)')
#plt.colorbar(term1,ax=axes.ravel().tolist(),ticks=[-10,0,10,20,30,40,50,60,70,80,90]) 
#plt.show()














