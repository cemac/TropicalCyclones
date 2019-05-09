########################################################################
# TC specific functions.                                               #
#                                                  John Ashcroft, 2017 #
########################################################################

import os
import matplotlib
#matplotlib.use('Agg')     #Use this when no graphic interface is available
import iris
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker





def find_centre_coords(md,TT,NN):
    filepath = "/nfs/a37/scjea/model_runs/Hagupit/u-af761/data/track_data/numpy_data/std1/track_data_{0:04d}{1:02d}00.npy.npz".format(md,TT)
    #print("=========== Data from u-af761 model run ===============")
    track_data = np.load(filepath)
    lats = track_data['arr_0']
    #print(lats)
    lons = track_data['arr_1']
    index = NN
    cenlat = lats[index]; cenlon = lons[index]
    return cenlat, cenlon

def find_centre_coords_numpy(md, TT, NN):
    filepath = "/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/{0:04d}_{1:02d}.npz".format(md,TT)
    track_data = np.load(filepath)
    lats = track_data['arr_0']
    lons = track_data['arr_1']
    index = NN / 6
    cenlat = lats[index]; cenlon = lons[index]
    return cenlat, cenlon


def find_centre_ens_4p4(md,TT,em,NN):
    filepath = "/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver_4p4_{0:04d}_{1:02d}Z_en{2:02d}.npz".format(md,TT,em)
    track_data = np.load(filepath)
    lats = track_data['lats']
    lons = track_data['lons']
    index = NN / 6
    cenlat = lats[index]; cenlon = lons[index]
    return cenlat, cenlon

def find_centre_ens_glm(md,TT,em,NN):
    filepath = "/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver_glm_{0:04d}_{1:02d}Z_en{2:02d}.npz".format(md,TT,em)
    track_data = np.load(filepath)
    lats = track_data['lats']
    lons = track_data['lons']
    index = NN / 6
    cenlat = lats[index]; cenlon = lons[index]
    return cenlat, cenlon

def find_centre_ens(md,TT,em,NN,mod):
    filepath = "/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver_{3}_{0:04d}_{1:02d}Z_en{2:02d}.npz".format(md,TT,em,mod)
    track_data = np.load(filepath)
    lats = track_data['lats']
    lons = track_data['lons']
    index = NN / 6
    cenlat = lats[index]; cenlon = lons[index]
    return cenlat, cenlon

def find_centre_3h(md,TT,em,v_day,v_time, mod):
    #if mod == '4p4':
    #    filepath = "/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver3hr_{3}_{0:04d}_{1:02d}Z_en{2:02d}.npz".format(md,TT,em,mod)
    #elif mod == 'glm':
    #
    df = "/nfs/a319/scjea/old_trackdata/orig_pctrack_firsttime_{0}_{1:02d}Z_{2}_em{3:02d}.npy".format(md,TT,mod,em)

#    if mod =='4p4':
#        filepath = "/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver3hr_{3}_{0:04d}_{1:02d}Z_en{2:02d}.npz".format(md,TT,em,mod)
#    elif mod == 'glm':
#        filepath = "/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver3hr_newinterp_{3}_{0:04d}_{1:02d}Z_en{2:02d}.npz".format(md,TT,em,mod)
    track_data = np.load(df).item()
    lats = track_data['lats']
    lons = track_data['lons']
    init_day = md % 100
    init_time = TT
    dt = (v_day - init_day)*24 + v_time - init_time
    index = dt/3
    #print(index)
#    days = track_data['vt_days']
#    times = track_data['vt_times']
    #index1 = np.where(days==v_day)
    #lats = lats[index1]; lons = lons[index1]
    #times = times[index1]
    #index2 = np.where(times==v_time)
    #cenlat = lats[index2]; cenlon = lons[index2]
    cenlat = lats[index]; cenlon = lons[index]
    return cenlat, cenlon


def load_tcverdata_3h_df(df0,storm,mod,md,TT):
    # Load all 12 ensemble members for the specified time and mdoel
    import matplotlib.dates as mdates
    import datetime as dt
    mth = int(str(md)[0:2]); day = int(str(md)[2:])
    ens_members = np.arange(45)
    i=0
    target_size = 39

    for em in ens_members:
        df = df0 + 'em{0:02d}.npz'.format(em)
        data = np.load(df)
        if storm == 'hagupit':
            #vt = np.arange(39)*3 + 18 + (md - 1201) * 24 + TT
            yr = 2014
        elif storm =='haiyan':
            yr = 2013
            #vt = np.arange(39)*3 + 18 + (md - 1103) * 24 + TT
        else:
            print('Storm unknown')
            exit()
        tcver_lats = data['lats']; tcver_lons = data['lons']; #tcver_press = data['minpress']; tcver_maxws = data['max_ws']
        tcver_days = data['vt_days']; tcver_hrs = data['vt_times']
        dates = []
        for day, hr in zip(tcver_days, tcver_hrs):
            date = dt.datetime(yr,mth,day,hr)
            date = mdates.date2num(date)
            dates.append(date)

        if i == 0:
            lats = tcver_lats; lons = tcver_lons; #press = tcver_press; maxws = tcver_maxws
            vts = dates
            i = i + 1
        else:
            lats = np.vstack([lats,tcver_lats]); lons = np.vstack([lons,tcver_lons]); #press = np.vstack([press,tcver_press]); maxws = np.vstack([maxws,tcver_maxws])
            vts = np.vstack([vts,dates])
    return lats,lons,vts

def track_data(filename):
    from netCDF4 import Dataset
    from netCDF4 import num2date
    dataset = Dataset(filename)
    lats =  dataset.variables['lat_wmo'][:]
    lons = dataset.variables['lon_wmo'][:]
    press = dataset.variables['pres_wmo'][:]
    mws = dataset.variables['wind_wmo'][:]
    times = num2date(dataset.variables['time_wmo'][:],dataset.variables['time_wmo'].units)
    times[:] = [time.replace(microsecond=0) for time in times]
    return lats,lons,press,mws,times
