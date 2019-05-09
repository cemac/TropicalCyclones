'''
Uses the data produced from my own tcver programme to plot various track and intensity
analysis results.
Requires improvement as for the CP where the storm leaves the domain (e.g. Haiyan and the
later forecasts) alot of the track errors diagnostics will fail.
'''


import os
#import matplotlib
#matplotlib.use('Agg')
import iris
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from netCDF4 import Dataset
from read_tctracker_data import *
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import itertools
from haversine_formula import haversine_d as calc_dist
from storm_tools import *
from cube_tools import *
from plotting_tools import *


def single_intensity(vts,fmw,storm,em):
    fmw = np.ma.masked_where(fmw==0,fmw)
    vts[:,1:] = np.ma.masked_where(vts[:,1:]==0,vts[:,1:])
    filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm)
    act_lats, act_lons, act_p, act_mws = track_data(filename)
    act_mws = act_mws * 0.51444444 #Knots to m/s
    time = np.arange(len(act_mws))*6
    colors = ['b','g','r','c','m','y','chocolate','b','g','r','c','m','y']
    linestyles = ['-','-','-','-','-','-','-','--','--','--','--','--','--']
    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    ax.plot(time,act_mws,'k-',linewidth=2,label='Best track data')
    plt.figure(1)
    plt.plot(vts[em,:],fmw[em,:],linewidth=1.5,linestyle='-',color='b')
    #### HAGUPIT AXES LABELLING
    if storm == 'Hagupit' or storm=='hagupit':
        intervals = np.arange(18,300,24)
        labels = ['01/12', '02/12', '03/12', '04/12', '05/12', '06/12', '07/12',  '08/12',  '09/12', '10/12', '11/12',  '12/12',  '13/12']
    #### HAIYAN AXES LABELLING
    elif storm == 'Haiyan':
        intervals = np.arange(18,300,24)
        labels = ['03/11', '04/11', '05/11', '06/11', '07/11', '07/11', '08/11',  '09/11',  '10/11', '11/11', '12/11',  '13/11',  '14/11']
    else:
        print 'Invalid storm name!'
        exit()
    ############################
    plt.xticks(intervals,labels,rotation='horizontal')
    plt.tick_params(labelsize=20)
    ax.set_xlabel('Date',fontsize='x-large')
    ax.set_ylabel('Windspeed, ms$^{-1}$', fontsize='x-large')
    plt.grid()
    plt.legend(loc='upper left')
    plt.show()
    plt.close()# Made better 8/3/2018

def single_track(lats,lons,storm,model,em):
    filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm)
    act_lats, act_lons, actp, actw = track_data(filename)
    Latitude=act_lats[:]; Longitude=act_lons[:]

    fig=plt.figure(figsize=(20,10), dpi=100)
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
    map_formatter(ax,labelsize=25,tick_base_x=5,tick_base_y=5)
    #ax.set_xlim([110,150])
    #ax.set_ylim([2.5, 16])
    ax.set_xlim([120,138])
    ax.set_ylim([7, 16])
    ######################### HAGUPIT LABELS  ############################
    if storm=='Hagupit' or storm=='hagupit':
        for i,j,k in zip(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],np.arange(len(Latitude[3:len(Latitude):4]))):
            if k==0:
                day=1
                month=12
                ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.4),xycoords='data', horizontalalignment='right',backgroundcolor='white',fontsize=20)
            else:
                day = k+1
                month=12
                if k<4: #4:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-1.0),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.2),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                #elif k==4:
                   # ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-1.0),xycoords='data', horizontalalignment='center',fontsize=20)
                #    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.8),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                #elif k==5:
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.5),xycoords='data', horizontalalignment='center',fontsize=20)
                #elif k==6:
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-1),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)#
                elif k==7:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i+0.2,j+0.6),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                else:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.6),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)#

    ############################### HAIYAN LABELS ##########################
    elif storm=='Haiyan' or storm == 'haiyan':
        for i,j,k in zip(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],np.arange(len(Latitude[3:len(Latitude):4]))):
            if k==0:
                day=3
                month=11
                ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.4),xycoords='data', horizontalalignment='right',backgroundcolor='white',fontsize=20)
            else:
                day = k+3
                month=11
                if k<5:
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.8),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-1.0),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                elif k==5:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.6),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                else:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.4),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)#
    #################################################################################

    colors = ['C0','C1','C2','C3','C4']
    labels = ['0-24h','24-48h','48-72h','72-96h','96-120h']
    lats = lats[em]; lons = lons[em]
    if lats.shape[0] == 39: # Missing T+0?
        for j in np.arange(5):
            if j==0:
                plt.plot(lons[0:8],lats[0:8],color=colors[j],linewidth=3,label=labels[j])
                
            elif j == 4:
                plt.plot(lons[j*8-1:],lats[j*8-1:],color=colors[j],linewidth=3,label=labels[j])
                
            else:
                plt.plot(lons[j*8-1:(j+1)*8],lats[j*8-1:(j+1)*8],color=colors[j],linewidth=3,label=labels[j])
    else:
        print "Error, check size of lat and lon datafiles"
        exit()
    best_track = ax.plot(Longitude,Latitude,'k-',linewidth=2.5,label='IBTrACS Best Track')
    daybesttrack = ax.plot(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],'gv',markersize=10,label='Position at 00Z')
    daybesttrack12 = ax.plot(Longitude[1:len(Longitude):4],Latitude[1:len(Latitude):4],'ko',markersize=6,label='Position at 12Z')

    plt.legend(loc=9,ncol=8,bbox_to_anchor=(0.5, -0.1))
    plt.show()
    plt.close()

def plot_track(lats,lons,plt_labels,storm,model): # Made better 8/3/18
    filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm)
    act_lats, act_lons, actp, actw = track_data(filename)
    Latitude=act_lats[:]; Longitude=act_lons[:]

    fig=plt.figure(figsize=(20,10), dpi=100)
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
    map_formatter(ax,labelsize=15,tick_base_x=10,tick_base_y=10)
    ax.set_xlim([110,155])
    ax.set_ylim([-2.5, 17.5])

    best_track = ax.plot(Longitude,Latitude,'k-',linewidth=2.5,label='IBTrACS Data')
    daybesttrack = ax.plot(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],'gv',markersize=10,label='Position at 0Z')

    ######################### HAGUPIT LABELS  ############################
    if storm=='Hagupit':
        for i,j,k in zip(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],np.arange(len(Latitude[3:len(Latitude):4]))):
            if k==0:
                day=1
                month=12
                ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.4),xycoords='data', horizontalalignment='right',size='x-large')
            else:
                day = k+1
                month=12
                if k<4:#k<5:
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.8),xycoords='data', horizontalalignment='center',size='x-large')
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.2),xycoords='data', horizontalalignment='center',size='x-large')
                elif k==5 or k==4:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.6),xycoords='data', horizontalalignment='center',size='x-large')
                else:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.4),xycoords='data', horizontalalignment='center',size='x-large')#

    ############################### HAIYAN LABELS ##########################
    elif storm=='Haiyan':
        for i,j,k in zip(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],np.arange(len(Latitude[3:len(Latitude):4]))):
            if k==0:
                day=3
                month=11
                ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.4),xycoords='data', horizontalalignment='right',size='x-large')
            else:
                day = k+3
                month=11
                if k<5:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.8),xycoords='data', horizontalalignment='center',size='x-large')
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.8),xycoords='data', horizontalalignment='center',size='x-large')
                elif k==5:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.6),xycoords='data', horizontalalignment='center',size='x-large')
                else:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.4),xycoords='data', horizontalalignment='center',size='x-large')#
    #################################################################################

    colors = ['b','g','r','c','m','y','chocolate','b','g','r','c','m','y']
    linestyles = ['-','-','-','-','-','-','-','--','--','--','--','--','--']

    for i in np.arange(lats.shape[0]):
        plt.plot(lons[i,:],lats[i,:],label=plt_labels[i],color=colors[i],linestyle=linestyles[i],linewidth=3)

    plt.legend(loc='lower left')
    plt.show()
    plt.close()

def plot_track_days(lats,lons,storm,model): # Made better 8/3/2018
    # TCVER produces data every 6 hours,therefore 4 connecting lines=24hrs
    # 5-day forecasts if the storm stays within the domain.
    filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm)
    act_lats, act_lons, actp, actw = track_data(filename)
    Latitude=act_lats[:]; Longitude=act_lons[:]

    fig=plt.figure(figsize=(20,10), dpi=100)
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
    map_formatter(ax,labelsize=25,tick_base_x=5,tick_base_y=5)
    ax.set_xlim([110,150])
    ax.set_ylim([2.5, 16])
    #ax.set_xlim([120,138])
    #ax.set_ylim([7, 16])


    ######################### HAGUPIT LABELS  ############################
    if storm=='Hagupit' or storm=='hagupit':
        for i,j,k in zip(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],np.arange(len(Latitude[3:len(Latitude):4]))):
            if k==0:
                day=1
                month=12
                ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.6),xycoords='data', horizontalalignment='right',backgroundcolor='white',fontsize=20,zorder=2)
            else:
                day = k+1
                month=12
                if k<4: #4:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.9),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20,zorder=2)
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.2),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                #elif k==4:
                   # ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-1.0),xycoords='data', horizontalalignment='center',fontsize=20)
                #    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.8),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                #elif k==5:
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.5),xycoords='data', horizontalalignment='center',fontsize=20)
                #elif k==6:
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-1),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)#
                elif k==7:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i+0.2,j+0.5),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20,zorder=2)
                else:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.5),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20,zorder=2)#

    ############################### HAIYAN LABELS ##########################
    elif storm=='Haiyan' or storm == 'haiyan':
        for i,j,k in zip(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],np.arange(len(Latitude[3:len(Latitude):4]))):
            if k==0:
                day=3
                month=11
                ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.4),xycoords='data', horizontalalignment='right',backgroundcolor='white',fontsize=20)
            else:
                day = k+3
                month=11
                if k<5:
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.8),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-1.3),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                elif k==5:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.8),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                else:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.6),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)#
    #################################################################################

    colors = ['C0','C1','C2','C3','C4']
    labels = ['0-24h','24-48h','48-72h','72-96h','96-120h']
    print lats.shape
    print lons.shape
    if lats.shape[1] == 20: # Missing T+0?
        for i in np.arange(lats.shape[0]):
            for j in np.arange(5):
                if j == 0:
                    if i == 0:
                        plt.plot(lons[i,0:4],lats[i,0:4],color=colors[j],linewidth=3,label='0-24h')
                    else:
                        plt.plot(lons[i,0:4],lats[i,0:4],color=colors[j],linewidth=3)
                else:
                    if i==0:
                        plt.plot(lons[i,j*4-1:(j+1)*4],lats[i,j*4-1:(j+1)*4],color=colors[j],linewidth=3,label=labels[j])
                    else:
                        plt.plot(lons[i,j*4-1:(j+1)*4],lats[i,j*4-1:(j+1)*4],color=colors[j],linewidth=3)
    elif lats.shape[1] == 21:
        for i in np.arange(lats.shape[0]):
            for j in np.arange(5):
                if i ==0:
                    plt.plot(lons[i,j*4:(j+1)*4],lats[i,j*4:(j+1)*4],color=colors[j],linewidth=3,label=labels[j])
                else:
                    plt.plot(lons[i,j*4:(j+1)*4],lats[i,j*4:(j+1)*4],color=colors[j],linewidth=3)
    elif lats.shape[1] == 39 or lats.shape[1] == 40: # Missing T+0 and T+120
        for i in np.arange(lats.shape[0]): #no. of ensemble members
            for j in np.arange(5):
                if j==0:
                    if i ==0:
                        plt.plot(lons[i,0:8],lats[i,0:8],color=colors[j],linewidth=3,label=labels[j],zorder=3)
                    else:
                        plt.plot(lons[i,0:8],lats[i,0:8],color=colors[j],linewidth=3,zorder=3)
                elif j == 4:
                    if i==0:
                        plt.plot(lons[i,j*8-1:],lats[i,j*8-1:],color=colors[j],linewidth=3,label=labels[j],zorder=3)
                    else:
                        plt.plot(lons[i,j*8-1:],lats[i,j*8-1:],color=colors[j],linewidth=3,zorder=3)
                else:
                    print j*8-1
                    if i==0:
                        plt.plot(lons[i,j*8-1:(j+1)*8],lats[i,j*8-1:(j+1)*8],color=colors[j],linewidth=3,label=labels[j],zorder=3)
                    else:
                        plt.plot(lons[i,j*8-1:(j+1)*8],lats[i,j*8-1:(j+1)*8],color=colors[j],linewidth=3,zorder=3)
    else:
        print "Error, check size of lat and lon datafiles"
        exit()

    best_track = ax.plot(Longitude,Latitude,'k-',linewidth=2.5,label='IBTrACS Best Track',zorder=1)
    daybesttrack = ax.plot(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],'v',color='k',markersize=20,label='Position at 00Z',zorder=4)
    daybesttrack12 = ax.plot(Longitude[1:len(Longitude):4],Latitude[1:len(Latitude):4],'ko',markersize=12,label='Position at 12Z',zorder=4)

    plt.legend(loc=9,ncol=8,bbox_to_anchor=(0.5, -0.1))
    plt.savefig('./{0}_{1}_track.png'.format(storm,model))
    plt.show()
    plt.close()

def calc_std_dev(lats,lons, md, TT,mod,include_storm=False):
    # Calculates the standard deviation of forecastes tracks, ignoring the storm
    # actual track. Feature to be added.
    t = []; std_dist = []; std_lat = []; std_lon = []
    num_times = len(lons[0,:])
    outfile='track_data/{0:04d}_{1:02d}Z_{2}_std'.format(md,TT,mod)
    for i in np.arange(num_times):
        dx = []
        t.append(i*6)
        avg_lon = np.mean(lons[:,i])
        avg_lat = np.mean(lats[:,i])
        for j in np.arange(len(lons[:,i])):
            dist = calc_dist(avg_lat,lats[j,i],avg_lon,lons[j,i])
            dx.append(dist)
        std_y = np.std(lats[:,i])
        std_x = np.std(lons[:,i])
        std_d = np.std(dx)
        std_dist.append(std_d); std_lat.append(std_y); std_lon.append(std_x)

    np.savez(outfile,t=t,std_dist=std_dist,std_lat=std_lat,std_lon=std_lon)    # Made better 8/3/2018

def plot_ws(vts,fmw,plt_labels,storm):
    fmw = np.ma.masked_where(fmw==0,fmw)
    vts[:,1:] = np.ma.masked_where(vts[:,1:]==0,vts[:,1:])
    filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm)
    act_lats, act_lons, act_p, act_mws = track_data(filename)
    act_mws = act_mws * 0.51444444 #Knots to m/s
    time = np.arange(len(act_mws))*6
    colors = ['b','g','r','c','m','y','chocolate','b','g','r','c','m','y']
    linestyles = ['-','-','-','-','-','-','-','--','--','--','--','--','--']

    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    ax.plot(time,act_mws,'k-',linewidth=2,label='Best track data')
    for i in np.arange(fmw.shape[0]):
        plt.figure(1)
        plt.plot(vts[i,:],fmw[i,:],label=plt_labels[i],linewidth=1.5)
    #### HAGUPIT AXES LABELLING
    if storm == 'Hagupit' or storm=='hagupit':
        intervals = np.arange(18,300,24)
        labels = ['01/12', '02/12', '03/12', '04/12', '05/12', '06/12', '07/12',  '08/12',  '09/12', '10/12', '11/12',  '12/12',  '13/12']
    #### HAIYAN AXES LABELLING
    elif storm == 'Haiyan':
        intervals = np.arange(18,300,24)
        labels = ['03/11', '04/11', '05/11', '06/11', '07/11', '07/11', '08/11',  '09/11',  '10/11', '11/11', '12/11',  '13/11',  '14/11']
    else:
        print 'Invalid storm name!'
        exit()
    ############################
    plt.xticks(intervals,labels,rotation='horizontal')
    plt.tick_params(labelsize=20)
    ax.set_xlabel('Date',fontsize='x-large')
    ax.set_ylabel('Windspeed, ms$^{-1}$', fontsize='x-large')
    plt.grid()
    plt.legend(loc='lower left',ncol=8)
    plt.show()
    plt.close()# Made better 8/3/2018

def plot_ws_box(vts, fmw, vtscp,fmwcp,storm):
    # Mask where the storm has left the domain
    # Requires both CP and global windspeed data
    fmw = np.ma.masked_where(fmw==0,fmw)
    fmwcp = np.ma.masked_where(fmwcp==0,fmwcp)
    vts[:,1:] = np.ma.masked_where(vts[:,1:]==0,vts[:,1:])
    vtscp[:,1:] = np.ma.masked_where(vtscp[:,1:]==0,vtscp[:,1:])
    avg_wind = []; ylow_err = []; yhigh_err = []
    avg_windcp = []; ylow_errcp = []; yhigh_errcp = []
    for i in np.arange(len(fmw[0,:])):
        avg_wind.append(np.mean(fmw[:,i]))
        ylow_err.append(np.mean(fmw[:,i]) - np.amin(fmw[:,i]))
        yhigh_err.append(np.amax(fmw[:,i]) - np.mean(fmw[:,i]))

    for i in np.arange(len(fmwcp[0,:])):
        avg_windcp.append(np.mean(fmwcp[:,i]))
        ylow_errcp.append(np.mean(fmwcp[:,i]) - np.amin(fmwcp[:,i]))
        yhigh_errcp.append(np.amax(fmwcp[:,i]) - np.mean(fmwcp[:,i]))

    fig = plt.figure(figsize=(15,10), dpi=100)
    ax = fig.add_subplot(1,1,1)
    filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm)
    act_lats, act_lons, act_p, act_mws = track_data(filename)
    act_mws = act_mws * 0.51444444
    time = np.arange(len(act_mws))*6
    ax.plot(time,act_mws,'k-',linewidth=2,label='Best track')

    plt.errorbar(vts[0,:],avg_wind,yerr=[ylow_err,yhigh_err],capsize=5,ecolor='C0')
    plt.plot(vts[0,:],avg_wind,color='C0',label='Global')
    plt.errorbar(vtscp[0,:],avg_windcp,yerr=[ylow_errcp,yhigh_errcp],capsize=5,color='C1',ecolor='C1')
    plt.plot(vtscp[0,:],avg_windcp,color='C1',label='Convection permitting')
    if storm == 'hagupit':
        #### HAGUPIT AXES LABELLING
        intervals = np.arange(18,300,24)
        labels = ['01/12', '02/12', '03/12', '04/12', '05/12', '06/12', '07/12',  '08/12',  '09/12', '10/12', '11/12',  '12/12',  '13/12']
    elif storm =='haiyan':
        #### HAIYAN AXES LABELLING
        intervals = np.arange(18,300,24)
        labels = ['03/11', '04/11', '05/11', '06/11', '07/11', '08/11', '09/11',  '10/11',  '11/11', '12/11', '13/11',  '14/11',  '15/11']
    else:
        print 'Invalid storm'
        exit()
    plt.xticks(intervals,labels,rotation='horizontal')
    plt.tick_params(labelsize=20)
    ax.set_xlabel('Date',fontsize=30)
    ax.set_ylabel('Windspeed, m$\,$s$^{-1}$', fontsize=30)
    ax.set_xlim([0,210])
    plt.grid()
    plt.legend(loc='upper left',fontsize=25)
    plt.show()

def plot_mslp(vts, fcp, plt_labels, storm):
    fcp = np.ma.masked_where(fcp==0,fcp)
    vts[:,1:] = np.ma.masked_where(vts[:,1:]==0,vts[:,1:])
    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm)
    act_lats, act_lons, act_p, act_mws = track_data(filename)
    time = np.arange(len(act_p))*6
    plt.plot(time,act_p,'k-',linewidth=2,label='Best track data')
    colors = ['b','g','r','c','m','y','chocolate','b','g','r','c','m','y']
    linestyles = ['-','-','-','-','-','-','-','--','--','--','--','--','--']
    if storm == 'Hagupit' or storm=='hagupit':
        intervals = np.arange(18,300,24)
        labels = ['01/12', '02/12', '03/12', '04/12', '05/12', '06/12', '07/12',  '08/12',  '09/12', '10/12', '11/12',  '12/12',  '13/12']
    elif storm == 'Haiyan':
        intervals = np.arange(18,300,24)
        labels = ['03/11', '04/11', '05/11', '06/11', '07/11', '08/11', '09/11', '10/11', '11/11']
    else:
        print 'Invalid storm'
        exit()
    print fcp.shape
    print vts.shape
    print len(plt_labels)
    print len(colors)
    for i in np.arange(fcp.shape[0]):
        plt.figure(1)
        plt.plot(vts[i,:],fcp[i,:],label=plt_labels[i],linewidth=1.5)
        #plt.savefig('{0}.png'.format(i))
        #plt.close()
    plt.xticks(intervals,labels,rotation='horizontal')
    plt.tick_params(labelsize=20)
    ax.set_xlabel('Date',fontsize='x-large')
    ax.set_ylabel('Minimum sea level pressure, hPa', fontsize='x-large')
    plt.grid()
    plt.legend(loc='lower left',ncol=8)
    plt.show()

def plot_mslp_box(vts, fmp, vtscp,fmpcp, storm):
    # Mask where the storm has left the domain
    # Requires both CP and global windspeed data
    fmp = np.ma.masked_where(fmp==0,fmp)
    fmpcp = np.ma.masked_where(fmpcp==0,fmpcp)
    vts[:,1:] = np.ma.masked_where(vts[:,1:]==0,vts[:,1:])
    vtscp[:,1:] = np.ma.masked_where(vtscp[:,1:]==0,vtscp[:,1:])
    avg_mslp = []; ylow_err = []; yhigh_err = []
    avg_mslpcp = []; ylow_errcp = []; yhigh_errcp = []
    for i in np.arange(len(fmp[0,:])):
        avg_mslp.append(np.mean(fmp[:,i]))
        ylow_err.append(np.mean(fmp[:,i]) - np.amin(fmp[:,i]))
        yhigh_err.append(np.amax(fmp[:,i]) - np.mean(fmp[:,i]))

    for i in np.arange(len(fmpcp[0,:])):
        avg_mslpcp.append(np.mean(fmpcp[:,i]))
        ylow_errcp.append(np.mean(fmpcp[:,i]) - np.amin(fmpcp[:,i]))
        yhigh_errcp.append(np.amax(fmpcp[:,i]) - np.mean(fmpcp[:,i]))

    fig = plt.figure(figsize=(15,10), dpi=100)
    ax = fig.add_subplot(1,1,1)
    filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm)
    act_lats, act_lons, act_p, act_mws = track_data(filename)
    time = np.arange(len(act_p))*6
    print time
    ax.plot(time,act_p,'k-',linewidth=2,label='Best track')

    plt.figure(1)
    plt.errorbar(vts[0,:],avg_mslp,yerr=[ylow_err,yhigh_err],capsize=5,ecolor='C0')
    plt.plot(vts[0,:],avg_mslp,color='C0',label='Global')
    plt.errorbar(vtscp[0,:],avg_mslpcp,yerr=[ylow_errcp,yhigh_errcp],capsize=5,color='C1',ecolor='C1')
    plt.plot(vtscp[0,:],avg_mslpcp,color='C1',label='Convection permitting')
    if storm == 'hagupit':
        #### HAGUPIT AXES LABELLING
        intervals = np.arange(18,300,24)
        labels = ['01/12', '02/12', '03/12', '04/12', '05/12', '06/12', '07/12',  '08/12',  '09/12', '10/12', '11/12',  '12/12',  '13/12']
    elif storm =='haiyan':
        #### HAIYAN AXES LABELLING
        intervals = np.arange(18,300,24)
        labels = ['03/11', '04/11', '05/11', '06/11', '07/11', '08/11', '09/11',  '10/11',  '11/11', '12/11', '13/11',  '14/11',  '15/11']
    else:
        print 'Invalid storm'
        exit()
    plt.xticks(intervals,labels,rotation='horizontal')
    plt.tick_params(labelsize=20)
    ax.set_xlabel('Date',fontsize=30)
    ax.set_ylabel('Sea level pressure (hPa)', fontsize=30)
    plt.grid()
    plt.legend(loc='lower left',fontsize=25)
    ax.set_xlim([0,210])
    plt.show()

def load_tcverdata(storm,mod,md,TT):
    # Load all 12 ensemble members for the specified time and mdoel
    ens_members = np.arange(12)
    i=0
    for em in ens_members:
        df = 'track_data/tcver3hr_new_{3}_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(md,TT,em,mod)
        data = np.load(df)
        if storm == 'hagupit':
            vt = np.arange(20)*6 + 18 + (md - 1201) * 24 + TT
        elif storm =='haiyan':
            vt = np.arange(20)*6 + 18 + (md - 1103) * 24 + TT
        else:
            print 'Storm unknown'
            exit()

        if i == 0:
            lats = data['lats']; lons = data['lons']; press = data['minpress']; maxws = data['max_ws']
            vts = vt
            i = i + 1
        else:
            lats = np.vstack([lats,data['lats']]); lons = np.vstack([lons,data['lons']]); press = np.vstack([press,data['minpress']]); maxws = np.vstack([maxws,data['max_ws']])
            vts = np.vstack([vts,vt])
    return lats,lons,press,maxws,vts



def load_tcverdata_3h_df(df0,storm,mod,md,TT):
    # Load all 12 ensemble members for the specified time and mdoel
    ens_members = np.arange(45)
    i=0
    target_size = 39
    for em in ens_members:
        df = df0 + 'em{0:02d}.npz'.format(em)
        data = np.load(df)
        if storm == 'hagupit' or storm=='Hagupit':
            vt = np.arange(39)*3 + 18 + (md - 1201) * 24 + TT
        elif storm =='haiyan':
            vt = np.arange(39)*3 + 18 + (md - 1103) * 24 + TT
        else:
            print 'Storm unknown'
            exit()
        tcver_lats = data['lats']; tcver_lons = data['lons']; tcver_press = data['minpress']; tcver_maxws = data['max_ws']
        print em
        while len(tcver_lats)<target_size:
            tcver_lats = np.append(tcver_lats, tcver_lats[-1])
            tcver_lons = np.append(tcver_lons, tcver_lons[-1])
            tcver_press = np.append(tcver_press,0)
            tcver_maxws = np.append(tcver_maxws,0)
        
        if i == 0:
            lats = tcver_lats; lons = tcver_lons; press = tcver_press; maxws = tcver_maxws
            vts = vt
            i = i + 1
        else:
            lats = np.vstack([lats,tcver_lats]); lons = np.vstack([lons,tcver_lons]); press = np.vstack([press,tcver_press]); maxws = np.vstack([maxws,tcver_maxws])
            vts = np.vstack([vts,vt])
    return lats,lons,press,maxws,vts


    

def load_tcverdata_3h_df_track(df0,storm,mod,md,TT):
    # Load all 12 ensemble members for the specified time and mdoel
    ens_members = np.arange(12)
    i=0
    target_size = 39
    for em in ens_members:
        df = df0 + 'em{0:02d}.npy'.format(em)
        data = np.load(df).item()
        if storm == 'hagupit':
            vt = np.arange(39)*3 + 18 + (md - 1201) * 24 + TT
        elif storm =='haiyan':
            vt = np.arange(39)*3 + 18 + (md - 1103) * 24 + TT
        else:
            print 'Storm unknown'
            exit()
        tcver_lats = data['lats']; tcver_lons = data['lons']; #tcver_press = data['minpress']; tcver_maxws = data['max_ws']
        while len(tcver_lats)<target_size:
            tcver_lats = np.append(tcver_lats, tcver_lats[-1])
            tcver_lons = np.append(tcver_lons, tcver_lons[-1])
            tcver_press = np.append(tcver_press,0)
            tcver_maxws = np.append(tcver_maxws,0)
        
        if i == 0:
            lats = tcver_lats; lons = tcver_lons; #press = tcver_press; maxws = tcver_maxws
            vts = vt
            i = i + 1
        else:
            lats = np.vstack([lats,tcver_lats]); lons = np.vstack([lons,tcver_lons]); #press = np.vstack([press,tcver_press]); maxws = np.vstack([maxws,tcver_maxws])
            vts = np.vstack([vts,vt])
    return lats,lons,vts


def load_tcverdata_3h(storm,mod,md,TT):
    # Load all 12 ensemble members for the specified time and mdoel
    ens_members = np.arange(12)
    i=0
    target_size = 39
    for em in ens_members:
        if mod == '4p4':
            df = 'track_data/tcver3hr_newinterp_{3}_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(md,TT,em,mod)
        else:
            df = 'track_data/tcver3hr_newinterp_{3}_{0:04d}_{1:02d}Z_en{2:02d}.npz'.format(md,TT,em,mod)
        data = np.load(df)
        if storm == 'hagupit':
            vt = np.arange(39)*3 + 18 + (md - 1201) * 24 + TT
        elif storm =='haiyan':
            vt = np.arange(39)*3 + 18 + (md - 1103) * 24 + TT
        else:
            print 'Storm unknown'
            exit()
        tcver_lats = data['lats']; tcver_lons = data['lons']; tcver_press = data['minpress']; tcver_maxws = data['max_ws']
        while len(tcver_lats)<target_size:
            tcver_lats = np.append(tcver_lats, tcver_lats[-1])
            tcver_lons = np.append(tcver_lons, tcver_lons[-1])
            tcver_press = np.append(tcver_press,0)
            tcver_maxws = np.append(tcver_maxws,0)
        
        if i == 0:
            lats = tcver_lats; lons = tcver_lons; press = tcver_press; maxws = tcver_maxws
            vts = vt
            i = i + 1
        else:
            lats = np.vstack([lats,tcver_lats]); lons = np.vstack([lons,tcver_lons]); press = np.vstack([press,tcver_press]); maxws = np.vstack([maxws,tcver_maxws])
            vts = np.vstack([vts,vt])
    return lats,lons,press,maxws,vts

def all_tracks(monthday, times, mod, storm):
    
    filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm)
    act_lats, act_lons, actp, actw = track_data(filename)
    Latitude=act_lats[:]; Longitude=act_lons[:]

    fig=plt.figure(figsize=(20,10), dpi=100)
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
    map_formatter(ax,labelsize=15,tick_base_x=10,tick_base_y=10)
    ax.set_xlim([100,155])
    ax.set_ylim([2.5, 25.])
    if storm == 'hagupit':
        labels = ['2/12, 00Z','2/12, 12Z','3/12, 00Z','3/12, 12Z','4/12, 00Z','4/12, 12Z','5/12, 00Z','5/12, 12Z','6/12, 00Z','6/12, 12Z','7/12, 00Z','7/12, 12Z']
    elif storm == 'haiyan':
        labels = ['4/11, 00Z','4/11, 12Z','5/11, 00Z','5/11, 12Z','6/11, 00Z','6/11, 12Z','7/11, 00Z','7/11, 12Z','8/11, 00Z','8/11, 12Z']
    
    J = 0
    for md in monthday:
        for TT in times:
            if md == 1207 and TT == 12:
                continue
            print [md,TT]
            [lats,lons,press,maxws,vts] = load_tcverdata_3h(storm,mod,md,TT)
            for i in np.arange(lats.shape[0]):
                for j in np.arange(1,len(lons[i,:])):
                    if lats[i,j] > 16 and lats[i,j-1] > lats[i,j] :
                        for k in np.arange(j,len(lons[i,:])):
                            lons[i,k] = lons[i,j-1]
                            lats[i,k] = lats[i,j-1]
                        break
                if J == 10:
                    if i == 0:
                        plt.plot(lons[i,:],lats[i,:],linewidth=1.5,color='chocolate',alpha=0.8,label=labels[J])
                    else:
                        plt.plot(lons[i,:],lats[i,:],linewidth=1.5,color='chocolate',alpha=0.8)
                else:
                    if i == 0:
                        plt.plot(lons[i,:],lats[i,:],linewidth=1.5,color='C{0}'.format(J),alpha=0.8,label=labels[J])
                    else:
                        plt.plot(lons[i,:],lats[i,:],linewidth=1.5,color='C{0}'.format(J),alpha=0.8)
            #avg_lat = np.mean(lats,axis=0)
            #avg_lon = np.mean(lons,axis=0)
            #plt.plot(avg_lon,avg_lat,linewidth=3,color='C{0}'.format(j))
            J = J+1
    best_track = ax.plot(Longitude,Latitude,'k-',linewidth=2.5)#,label='IBTrACS Data'
    daybesttrack = ax.plot(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],'gv',markersize=10) #,label='Position at 0Z'
    plt.legend(loc=9,ncol=11,bbox_to_anchor=(0.5, -0.1))
    plt.show()
        
def plot_errors(md,TT,mod,ens_members):
    fig = plt.figure(figsize=(20,10),dpi=100)
    ax = fig.add_subplot(111)
    times = np.arange(19) * 6 + 6
    colors = ['C0','C1','C2','C3','C4','C5','C0','C1','C2','C3','C4','C5']
    linestyles = ['-','-','-','-','-','-','--','--','--','--','--','--']
    for em in ens_members:
        df = 'track_data/trackerrors_{0}_{1:02d}_{2}_{3:02d}.npz'.format(md,TT,mod,em)
        data = np.load(df)
        dpe = data['dpe']/1000.; cte = data['cte']/1000.; ate = data['ate']/1000.
        ax.plot(times,dpe,linewidth=2,linestyle = linestyles[em], color=colors[em],label='Ensemble member {0:02d}'.format(em))
        
#        ax.plot(np.absolute(cte))
#        ax.plot(ate)
    ax.set_ylabel('DPE (km)',fontsize=30)
    ax.set_xlabel('Forecast time (hrs)',fontsize=30)
    plt.grid()
    plt.legend(fontsize=15)
    plt.tick_params(labelsize=20)
    plt.show()
    
    plt.close()
    
def avg_errors():
    fig = plt.figure(figsize=(20,10),dpi=100)
    ax = fig.add_subplot(111)
    times = np.arange(19)*6 + 6
    # Hag, glm:
    ens_members = np.arange(12)
    for em in ens_members:
        df = 'track_data/trackerrors_1203_12_glm_{0:02d}.npz'.format(em)
        data = np.load(df)
        dpe = data['dpe']/1000.; cte = data['cte']/1000.; ate = data['ate']/1000.
        if em == 0:
            dpe_list = dpe; ate_list = ate; cte_list = cte
        else:
            dpe_list = np.vstack([dpe_list,dpe]); ate_list = np.vstack([ate_list,ate])
            cte_list = np.vstack([cte_list,cte])
        
    hag_glm_dpe = np.mean(dpe_list,axis=0); hag_glm_ate = np.mean(np.absolute(ate_list),axis=0)
    hag_glm_cte = np.mean(np.absolute(cte_list),axis=0)
    
    for em in ens_members:
        df = 'track_data/trackerrors_1203_12_4p4_{0:02d}.npz'.format(em)
        data = np.load(df)
        dpe = data['dpe']/1000.; cte = data['cte']/1000.; ate = data['ate']/1000.
        if em == 0:
            dpe_list = dpe; ate_list = ate; cte_list = cte
        else:
            dpe_list = np.vstack([dpe_list,dpe]); ate_list = np.vstack([ate_list,ate])
            cte_list = np.vstack([cte_list,cte])
        
    hag_4p4_dpe = np.mean(dpe_list,axis=0); hag_4p4_ate = np.mean(np.absolute(ate_list),axis=0)
    hag_4p4_cte = np.mean(np.absolute(cte_list),axis=0)
    
    for em in ens_members:
        df = 'track_data/trackerrors_1104_12_glm_{0:02d}.npz'.format(em)
        data = np.load(df)
        dpe = data['dpe']/1000.; cte = data['cte']/1000.; ate = data['ate']/1000.
        if em == 0:
            dpe_list = dpe; ate_list = ate; cte_list = cte
        else:
            dpe_list = np.vstack([dpe_list,dpe]); ate_list = np.vstack([ate_list,ate])
            cte_list = np.vstack([cte_list,cte])
        
    hai_glm_dpe = np.mean(dpe_list,axis=0); hai_glm_ate = np.mean(np.absolute(ate_list),axis=0)
    hai_glm_cte = np.mean(np.absolute(cte_list),axis=0)

    for em in ens_members:
        df = 'track_data/trackerrors_1104_12_4p4_{0:02d}.npz'.format(em)
        data = np.load(df)
        dpe = data['dpe']/1000.; cte = data['cte']/1000.; ate = data['ate']/1000.
        if em == 0:
            dpe_list = dpe; ate_list = ate; cte_list = cte
        else:
            dpe_list = np.vstack([dpe_list,dpe]); ate_list = np.vstack([ate_list,ate])
            cte_list = np.vstack([cte_list,cte])
        
    hai_4p4_dpe = np.mean(dpe_list,axis=0); hai_4p4_ate = np.mean(np.absolute(ate_list),axis=0)
    hai_4p4_cte = np.mean(np.absolute(cte_list),axis=0)
    
    ax.plot(times,hag_glm_dpe,linestyle='-',color='C0',label='Hagupit, global ensemble',linewidth=2)
    ax.plot(times,np.absolute(hag_glm_ate),linestyle='--',color='C0',linewidth=2)
    ax.plot(times,np.absolute(hag_glm_cte),linestyle=':',color='C0',linewidth=2)
    ax.plot(times,hag_4p4_dpe,linestyle='-',color='C1',label='Hagupit, CP ensemble',linewidth=2)
    ax.plot(times,np.absolute(hag_4p4_ate),linestyle='--',color='C1',linewidth=2)
    ax.plot(times,np.absolute(hag_4p4_cte),linestyle=':',color='C1',linewidth=2)
    ax.plot(times,hai_glm_dpe,linestyle='-',color='C2',label='Haiyan, global ensemble',linewidth=2)
    ax.plot(times,np.absolute(hai_glm_ate),linestyle='--',color='C2',linewidth=2)
    ax.plot(times,np.absolute(hai_glm_cte),linestyle=':',color='C2',linewidth=2)
    ax.plot(times,hai_4p4_dpe,linestyle='-',color='C3',label='Haiyan, CP ensemble',linewidth=2)
    ax.plot(times,np.absolute(hai_4p4_ate),linestyle='--',color='C3',linewidth=2)
    ax.plot(times,np.absolute(hai_4p4_cte),linestyle=':',color='C3',linewidth=2)
    
    ax.set_ylabel('Error (km)',fontsize=30)
    ax.set_xlabel('Forecast time (hrs)',fontsize=30)
    plt.grid()
    plt.legend(fontsize=15)
    plt.tick_params(labelsize=20)
    plt.show()
    
    
def plot_mslp_box_singleforecast(vts, fmp, ax,label,color_code):
    # Mask where the storm has left the domain
    # Requires ax to be added as input
    fmp = np.ma.masked_where(fmp==0,fmp)
    vts[:,1:] = np.ma.masked_where(vts[:,1:]==0,vts[:,1:])
    avg_mslp = []; ylow_err = []; yhigh_err = []
    for i in np.arange(len(fmp[0,:])):
        avg_mslp.append(np.mean(fmp[:,i]))
        ylow_err.append(np.mean(fmp[:,i]) - np.amin(fmp[:,i]))
        yhigh_err.append(np.amax(fmp[:,i]) - np.mean(fmp[:,i]))

    ax.errorbar(vts[0,:],avg_mslp,yerr=[ylow_err,yhigh_err],capsize=5,ecolor=color_code, color=color_code)
    ax.plot(vts[0,:],avg_mslp,color=color_code,label=label)

def fortis_track_days(lats,lons,storm,model,act_lats,act_lons): # Made better 8/3/2018
    # TCVER produces data every 6 hours,therefore 4 connecting lines=24hrs
    # 5-day forecasts if the storm stays within the domain.
    filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm)
    act_lats, act_lons, actp, actw = track_data(filename)
    Latitude=act_lats[:]; Longitude=act_lons[:]

    fig=plt.figure(figsize=(20,10), dpi=100)
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
    map_formatter(ax,labelsize=25,tick_base_x=5,tick_base_y=5)
    ax.set_xlim([120,140])
    ax.set_ylim([5, 20])
    #ax.set_xlim([110,140])
    #ax.set_ylim([5.875, 16])
    #ax.set_xlim([120,138])
    #ax.set_ylim([7, 16])


    
    colors = ['C0','C1','C2','C3','C4']
    labels = ['0-24h','24-48h','48-72h','72-96h','96-120h']
    sinlge_lat = lats[0]; single_lons = lats[0]
    if lats.shape[1] == 20: # Missing T+0?
        for j in np.arange(5):
            if j == 0:
                plt.plot(lons[0,0:4],lats[0,0:4],color=colors[j],linewidth=3)
            else:
                plt.plot(lons[0,j*4-1:(j+1)*4],lats[0,j*4-1:(j+1)*4],color=colors[j],linewidth=3)
        for i in np.arange(lats.shape[0]):
            for j in np.arange(5):
                if j == 0:
                    if i == 0:
                        plt.plot(lons[i,0:4],lats[i,0:4],color=colors[j],linewidth=3,label='0-24h')
                    else:
                        plt.plot(lons[i,0:4],lats[i,0:4],color=colors[j],linewidth=3)
                else:
                    if i==0:
                        plt.plot(lons[i,j*4-1:(j+1)*4],lats[i,j*4-1:(j+1)*4],color=colors[j],linewidth=3,label=labels[j])
                    else:
                        plt.plot(lons[i,j*4-1:(j+1)*4],lats[i,j*4-1:(j+1)*4],color=colors[j],linewidth=3)
    elif lats.shape[1] == 21:
        for i in np.arange(lats.shape[0]):
            for j in np.arange(5):
                if i ==0:
                    plt.plot(lons[i,j*4:(j+1)*4],lats[i,j*4:(j+1)*4],color=colors[j],linewidth=3,label=labels[j])
                else:
                    plt.plot(lons[i,j*4:(j+1)*4],lats[i,j*4:(j+1)*4],color=colors[j],linewidth=3)
    elif lats.shape[1] == 39: # Missing T+0 and T+120
        for j in np.arange(5):
            if j==0:
                plt.plot(lons[0,0:8],lats[0,0:8],color=colors[j],linewidth=3)
            elif j == 4:
                plt.plot(lons[0,j*8-1:],lats[0,j*8-1:],color=colors[j],linewidth=3)
            else:
                plt.plot(lons[0,j*8-1:(j+1)*8],lats[0,j*8-1:(j+1)*8],color=colors[j],linewidth=3)
        #plt.savefig('./{0}_{1}_tracks1.png'.format(storm,model))
        for i in np.arange(lats.shape[0]): #no. of ensemble members
            for j in np.arange(5):
                if j==0:
                    if i ==0:
                        plt.plot(lons[i,0:8],lats[i,0:8],color=colors[j],linewidth=3,label=labels[j])
                    else:
                        plt.plot(lons[i,0:8],lats[i,0:8],color=colors[j],linewidth=3)
                elif j == 4:
                    if i==0:
                        plt.plot(lons[i,j*8-1:],lats[i,j*8-1:],color=colors[j],linewidth=3,label=labels[j])
                    else:
                        plt.plot(lons[i,j*8-1:],lats[i,j*8-1:],color=colors[j],linewidth=3)
                else:
                    print j*8-1
                    if i==0:
                        plt.plot(lons[i,j*8-1:(j+1)*8],lats[i,j*8-1:(j+1)*8],color=colors[j],linewidth=3,label=labels[j])
                    else:
                        plt.plot(lons[i,j*8-1:(j+1)*8],lats[i,j*8-1:(j+1)*8],color=colors[j],linewidth=3)
    elif lats.shape[1] == 40: # Missing T+0 and T+120
        for i in np.arange(lats.shape[0]): #no. of ensemble members
            for j in np.arange(5):
                if j==0:
                    if i ==0:
                        plt.plot(lons[i,0:9],lats[i,0:9],color=colors[j],linewidth=3,label=labels[j])
                    else:
                        plt.plot(lons[i,0:9],lats[i,0:9],color=colors[j],linewidth=3)
                elif j == 4:
                    if i==0:
                        plt.plot(lons[i,j*8:],lats[i,j*8:],color=colors[j],linewidth=3,label=labels[j])
                    else:
                        plt.plot(lons[i,j*8:],lats[i,j*8:],color=colors[j],linewidth=3)
                else:
                    if i==0:
                        plt.plot(lons[i,j*8:(j+1)*8+1],lats[i,j*8:(j+1)*8+1],color=colors[j],linewidth=3,label=labels[j])
                    else:
                        plt.plot(lons[i,j*8:(j+1)*8+1],lats[i,j*8:(j+1)*8+1],color=colors[j],linewidth=3)
        
    else:
        print "Error, check size of lat and lon datafiles"
        print lats.shape
        exit()

    best_track = ax.plot(Longitude,Latitude,'k-',linewidth=3.5,label='IBTrACS Best Track')
    daybesttrack = ax.plot(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],'gv',markersize=10,label='Position at 00Z')
    daybesttrack12 = ax.plot(Longitude[1:len(Longitude):4],Latitude[1:len(Latitude):4],'ko',markersize=8,label='Position at 12Z')
    ######################### HAGUPIT LABELS  ############################
    if storm=='Hagupit' or storm=='hagupit':
        for i,j,k in zip(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],np.arange(len(Latitude[3:len(Latitude):4]))):
            if k==0:
                day=1
                month=12
                ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.4),xycoords='data', horizontalalignment='right',backgroundcolor='white',fontsize=20)
            else:
                day = k+1
                month=12
                if k<4: #4:
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-1.0),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.4),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                #elif k==4:
                   # ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-1.0),xycoords='data', horizontalalignment='center',fontsize=20)
                #    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.8),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                #elif k==5:
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.5),xycoords='data', horizontalalignment='center',fontsize=20)
                #elif k==6:
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-1),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)#
                elif k==7:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i+0.2,j+0.4),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                else:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.4),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)#

    ############################### HAIYAN LABELS ##########################
    elif storm=='Haiyan' or storm == 'haiyan':
        for i,j,k in zip(Longitude[3:len(Longitude):4],Latitude[3:len(Latitude):4],np.arange(len(Latitude[3:len(Latitude):4]))):
            if k==0:
                day=3
                month=11
                ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.4),xycoords='data', horizontalalignment='right',backgroundcolor='white',fontsize=20)
            else:
                day = k+3
                month=11
                if k<5:
                    #ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-0.8),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j-1.0),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                elif k==5:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.6),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)
                else:
                    ax.annotate('{0:02d}/{1:02d}'.format(day,month),xy=(i,j+0.4),xycoords='data', horizontalalignment='center',backgroundcolor='white',fontsize=20)#
    #################################################################################


    plt.legend(loc=9,ncol=8,bbox_to_anchor=(0.5, -0.1))
    plt.show()
    #plt.savefig('./{0}_{1}_tracks3.png'.format(storm,model))
    plt.close()

     
def main():
################################################################################
    monthday = [1104,1105,1106,1107,1108]
    #monthday = [1203,1204,1205,1206]
    #monthday = [1203]; times = [12]
    monthday = [1104]
    times = [12]
    mod = 'glm'; storm = 'haiyan'
    labels = ['Ensemble member {0:02d}'.format(i) for i in np.arange(45)]
    #avg_errors()
    for md in monthday:
        for TT in times:
            filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm)
            act_lats, act_lons, act_p, act_mws = track_data(filename)
            #df = 'track_data/tcver3hr_{2}_{0:04d}_{1:02d}Z_'.format(md,TT,mod)
            df = '/nfs/a319/scjea/old_trackdata/orig_pctrack_{0}_{1:02d}Z_{2}_'.format(md,TT,mod)
            [lats,lons,vts] = load_tcverdata_3h_df_track(df,storm,mod,md,TT)
            plot_track_days(lats,lons,storm,mod)
            #mod = '4p4'
            #df = 'track_data/pctracker_{2}_{0:04d}_{1:02d}Z_'.format(md,TT,mod)
            #[latscp,lonscp,presscp,maxwscp,vts] = load_tcverdata_3h_df(df,storm,mod,md,TT)
            #plot_ws_box(vts,maxws,vts,maxwscp,'hagupit')
            #plot_ws(vts,maxws,labels,'hagupit')
            #df = 'track_data/pctracker_bigens_{0}_{1}_{2:02d}_'.format(mod,md,TT)
            #[lats,lons,vts] = load_tcverdata_3h_df_track(df,storm,mod,md,TT)
            #print vts
            #fortis_track_days(lats,lons,storm,mod,act_lats,act_lons)
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            #print 'Loading data for {0}, md={1}, TT={2:02d}'.format(storm,md,TT)
            #[lats,lons,press,maxws,vts] = load_tcverdata_3h(storm,'glm',md,TT)
            ##plot_track_days(lats1,lons1,storm,mod)
            #[latscp,lonscp,presscp,maxwscp,vtscp] = load_tcverdata_3h(storm,'4p4',md,TT)
            ##print lats[0].shape[0]
         ##   single_track(lats,lons,storm,mod,em)
            ##single_intensity(vts,maxws,storm,em)
            ##[lats,lons,presscp,maxwscp,vtscp] = load_tcverdata(storm,'4p4',md,TT)
##            indexs = np.arange(20) * 2 
##            press = press[:,indexs]; presscp = presscp[:,indexs]; vts=vts[:,indexs]; vtscp=vtscp[:,indexs]
##            maxws = maxws[:,indexs]; maxwscp = maxwscp[:,indexs]
            
            ##plot_ws_box(vts, maxws, vtscp,maxwscp,storm)
            #plot_mslp_box(vts,press,vtscp,presscp,storm)
            
            #plot_track_days(lats,lons,storm,mod)
            #print lats1[0]
            #print lats[0]
            #plot_track_days(lats,lons,storm,mod)
            #calc_std_dev(lats,lons,md,TT,mod)












########################### Still to edit: #####################################


def plot_poserror_dpe1(ft,dpe,plt_labels):
    dpe = np.ma.masked_where(dpe==0,dpe)
    ft[:,1:-1] = np.ma.masked_where(ft[:,1:-1]==0,ft[:,1:-1])
    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    colors = ['b','g','r','c','m','y','chocolate','b','g','r','c','m','y']
    linestyles = ['-','-','-','-','-','-','-','--','--','--','--','--','--']

    for i in np.arange(dpe.shape[0]):
        plt.figure(1)
        vt = ft[i,:]
        plt.plot(vt,dpe[i,:],label=plt_labels[i],linewidth=2,linestyle=linestyles[i],color=colors[i])
    plt.grid()
    ax.set_xlabel('Forecast time (hr)',fontsize='x-large')
    ax.set_ylabel('Direct poistional error, km', fontsize='x-large')
    plt.legend(loc='upper left')
    plt.show()


def plot_poserror_dpe2(ft,dpe,plt_labels):
    dpe = np.ma.masked_where(dpe==0,dpe)
    ft[:,1:-1] = np.ma.masked_where(ft[:,1:-1]==0,ft[:,1:-1])
    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    colors = ['b','g','r','c','m','y','chocolate','b','g','r','c','m','y']
    linestyles = ['-','-','-','-','-','-','-','--','--','--','--','--','--']

    for i in np.arange(dpe.shape[0]):
        plt.figure(1)
        vt = ft[i,:] + (i+1)*12
        plt.plot(vt,dpe[i,:],label=plt_labels[i],linewidth=2,color=colors[i],linestyle=linestyles[i])
    plt.xlim(xmin=0)
    intervals = np.arange(12,300,24)
    labels = ['01/12', '02/12', '03/12', '04/12', '05/12', '06/12', '07/12',  '08/12',  '09/12', '10/12', '11/12',  '12/12',  '13/12']
    plt.xticks(intervals,labels,rotation='horizontal')
    plt.tick_params(labelsize=20)
    ax.set_xlabel('Date',fontsize='x-large')
    ax.set_ylabel('Direct positional error, km', fontsize='x-large')
    plt.grid()
    plt.legend(loc='upper left')
    plt.show()

def make_same_size(array, sarray):
    if array.shape == sarray.shape:
        return sarray
    else:
        sarray = np.append(sarray,0)
        return sarray

def crop_data(data):
    # Lat/lons may have been computed to keep the last known position as the storm goes off the domain. This function changes these positions to 0 ahead of being masked in the track function
    ## Haiyan domain --> (()) (())
    for i in np.arange(len(data) - 1):
        if data[i] == data[-1]:
            for j in np.arange(i+1, len(data), 1):
                data[j] = 0
    return data

def track_data(filename):
    dataset = Dataset(filename)
    print dataset
    lats =  dataset.variables['lat_for_mapping'][:]
    lons = dataset.variables['lon_for_mapping'][:]
    press = dataset.variables['pres_wmo'][:]
    mws = dataset.variables['wind_wmo'][:]
    return lats,lons,press,mws








if __name__ == '__main__':
    main()
