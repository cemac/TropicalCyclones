#######################################################################
# TC tracker using UM data. Technique based on that used in the TCVER #
# program by Julian Hemming.                                          #
#                                            John Ashcroft, Sept 2017 #
#######################################################################

# Edit - May, 2018
##  Added interpolation onto finer grid for grid spaces > 0.08 (chosen arbritarily)

import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from netCDF4 import Dataset
from cube_tools import *
import iris
import iris.analysis.calculus
import iris.analysis
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from interp_to_grid import *

def main():
############################## Change this #####################################
    monthday = [1105,1106,1107,1108]; times = [12,0]
    ens_members = np.arange(12)
#    ens_members = [4]
    ens_members = np.arange(12)
    models = ['4p4']
    storm = 'Haiyan'
    month = 11
    year = 2013
    dt = 3 #hours
################################################################################
    if storm == 'Hagupit':
        storm_id = 'u-ap087'
    elif storm=='Haiyan':
        storm_id = 'u-ap121'

    for md in monthday:
        for TT in times:
            for mod in models:
                if storm == 'Hagupit':
                    init_day = md - 1200
                elif storm == 'Haiyan':
                    init_day = md - 1100
                init_time = TT
                print '***** Finding initial search position. ******'
                lat0, lon0 = find_init_sc2(md, TT, storm)
                lat1, lon1 = find_init_sc2(md, TT+6,storm)
                #lat0 is at time t+0, lat1 at time t+6. We want time t+3 and t+6 so....
                lat0 = 0.5*(lat1+lat0)
                lon0 = 0.5*(lon1+lon0)
                print 'Initial search position: ({0:.2f},{1:.2f})'.format(lat0,lon0)
                print '***** Beginning main tcver program. *****'
                for em in ens_members:
                    df = '/nfs/a37/scjea/model_runs/{0}/{1}/data/{2}/{3}_{4:02d}Z/tctracker/tcver3hr_{2}_{3}_{4:02d}Z_em{5:02d}.pp'.format(storm,storm_id,mod,md,TT,em)
                    print 'Ensemble member: {0:02d}'.format(em)
                    lats, lons, minpress, max_ws,valid_days, valid_times = tcver(df,init_day,init_time,lat0,lon0,lat1,lon1,yr=year,mth=month)
                #while len(lats)<20:
            #        lats.append(None); lons.append(None); minpress.append(None); max_ws.append(None)
                    np.savez('./track_data/tcver3hr_newinterp_{3}_{0:04d}_{1:02d}Z_en{2:02d}'.format(md,TT,em,mod),lats=lats,lons=lons,minpress=minpress,max_ws=max_ws,vt_days=valid_days,vt_times=valid_times)
            print '***** End of program *****'

def closed_isobar_check(slp,cenlat,cenlon,minslp,rad=3):
    phis = np.arange(0,2*np.pi,np.pi / 4.)
    storm_centre = True
    for phi in phis:
        lonpoi = cenlon + rad * np.cos(phi)
        latpoi = cenlat + rad * np.sin(phi)
        poi = [('longitude', lonpoi), ('latitude',latpoi)]
        slp_check = slp.interpolate(poi,iris.analysis.Linear())
        if slp_check.data < minslp + 100.0:
            storm_centre = False
    return storm_centre


def tcver(df,init_day,init_time,lat0,lon0,lat1,lon1,yr=2014,mth=12):
    # 5 day forecasts - so we wun from init_day to init day + 5. However ignore the last time
    dt = 3
    days = np.arange(init_day, init_day + 6, 1)
    times = np.arange(8)*3  # i.e. every 3 hours.
    # Define the search radius for vorticity and slp - for 3 hour forecast use 1.5
    vrtrad = 3.
    slprad = 3.
    # Stash codes to load variables
    u_stash = 'm01s15i201'
    v_stash = 'm01s15i202'
    sfcu_stash = 'm01s03i225'
    sfcv_stash = 'm01s03i226'
    slp_stash = 'm01s16i222'

    plev = 850 # Data should only be on 850 hpa anyway, but just in case...
    # Constraints for cubes:
    pressure_constraint = iris.Constraint(pressure=plev)
    u_constraint = iris.AttributeConstraint(STASH=u_stash)
    v_constraint = iris.AttributeConstraint(STASH=v_stash)
    slp_constraint = iris.AttributeConstraint(STASH=slp_stash)
    sfcu_constraint = iris.AttributeConstraint(STASH=sfcu_stash)
    sfcv_constraint = iris.AttributeConstraint(STASH=sfcv_stash)

    # Load the cubes at all times
    uwind = iris.load_cube(df,u_constraint).extract(pressure_constraint)
    vwind = iris.load_cube(df,v_constraint).extract(pressure_constraint)
    mslp = iris.load_cube(df,slp_constraint)
    sfc_uwind = iris.load_cube(df,sfcu_constraint)
    sfc_vwind = iris.load_cube(df,sfcv_constraint)

    # Find coordintes of domain
    maxlat = np.amax(uwind.coord('latitude').points)
    minlat = np.amin(uwind.coord('latitude').points)
    maxlon = np.amax(uwind.coord('longitude').points)
    minlon = np.amin(uwind.coord('longitude').points)
    # Initialise arrays to populate
    latitudes = []; longitudes = []; minpress = []; maxws = []
    valid_days = []; valid_times = []
    cont_tracking = True
    for dd in days:
        for hr in times:
            # Loop through every time, ignoring invalid times.
            if dd < init_day or (dd == init_day and hr < init_time + 3): # may have to also ignore T+0 <<-- we do as no slp
                continue
            if dd > init_day + 5 or (dd == init_day + 5 and (hr > init_time - dt or init_time == 0)):
                continue
            # Load data at the current time
            time_constraint = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr))
            with iris.FUTURE.context(cell_datetime_objects=True):
                u = uwind.extract(time_constraint)
                v = vwind.extract(time_constraint)
                sfc_u = sfc_uwind.extract(time_constraint)
                sfc_v = sfc_vwind.extract(time_constraint)
                slp = mslp.extract(time_constraint)
            # 10m windspeed...
            ws = (sfc_u**2 + sfc_v**2)**0.5
            # 850 hPa vorticity...
            vrt = iris.analysis.calculus.curl(u,v)[2]
            if dd == init_day and hr == init_time+3: # i.e. first timesteo available, initial search position uses the observed position
                b_constraint = box_constraint(lat0-vrtrad,lat0+vrtrad,lon0-vrtrad,lon0+vrtrad)
            elif dd == init_day and hr == init_time + 6: # i.e. TS2, search position still uses the observed location
                b_constraint = box_constraint(lat1-vrtrad,lat1+vrtrad,lon1-vrtrad,lon1+vrtrad)
            else:
                lat2, lon2 = interp_coords(latitudes[-2],longitudes[-2],latitudes[-1],longitudes[-1]) # Interpolation to define next search position
                if lat2 + vrtrad > maxlat or lat2 - vrtrad < minlat or lon2 + vrtrad > maxlon or lon2 - vrtrad< minlon:
                    print 'WARNING: Interpolated position of storm search box out of domain, exiting tcver'
                    return latitudes, longitudes, minpress, maxws, valid_days, valid_times
                b_constraint = box_constraint(lat2-vrtrad,lat2+vrtrad,lon2-vrtrad,lon2+vrtrad)
            # Reduce the cube to save time
            reduced_vrt = vrt.extract(b_constraint)
            # Find the max vorticity
            maxvrt_lat, maxvrt_lon, maxvrt = max_vals(reduced_vrt, rad=vrtrad)
            # define new position according to the max vorticity gridpoint
            if maxvrt_lat + slprad > maxlat or maxvrt_lat - slprad < minlat or maxvrt_lon + slprad > maxlon or maxvrt_lon - slprad< minlon:
                print 'WARNING: Interpolated position of storm search box out of domain, exiting tcver'
                return latitudes, longitudes, minpress, maxws, valid_days, valid_times
            b_constraint = box_constraint(maxvrt_lat-slprad, maxvrt_lat+slprad, maxvrt_lon-slprad, maxvrt_lon+slprad)
            reduced_slp = slp.extract(b_constraint)
            # Find the mslp gridpoint
            minslp_lat, minslp_lon, min_slp = min_vals(reduced_slp, rad=slprad)
            if cont_tracking == True:
                cont_tracking = closed_isobar_check(reduced_slp,minslp_lat, minslp_lon, min_slp)
            if cont_tracking == False:
                print 'WARNING: Closed isobar check failed.'
                return latitudes, longitudes, minpress, maxws, valid_days, valid_times
            #print min_slp
            b_constraint = box_constraint(minslp_lat-slprad,minslp_lat+slprad,minslp_lon-slprad,minslp_lon+slprad)
            # Reduce the windspeed to just near the mslp gridpoint
            reduced_ws = ws.extract(b_constraint)
            # Find dlon and dlat
            dlon = abs(reduced_slp.coord('longitude').points[5] - reduced_slp.coord('longitude').points[4])
            dlat = abs(reduced_slp.coord('latitude').points[5] - reduced_slp.coord('latitude').points[4])
            if dlon > 0.00008: # i.e. coarse global grid # everything for now
                slprad2a = dlat*2.5; slprad2b = dlon*2.5
                if minslp_lat + slprad2a > maxlat or minslp_lat - slprad2a < minlat or minslp_lon + slprad2b > maxlon or minslp_lon - slprad2b< minlon:
                    print 'WARNING: Interpolated position of storm search box out of domain, exiting tcver'
                    return latitudes, longitudes, minpress, maxws, valid_days, valid_times
                b_constraint = box_constraint(minslp_lat-2.5*dlat, minslp_lat+2.5*dlat,minslp_lon-2.5*dlon,minslp_lon+2.5*dlon)
                reduced_slp2 = reduced_slp.extract(b_constraint)
                # new slp reduction ready for interpolation
                # Initial technique: ## Disused in favour of next technique - see notes, John Ashcroft, June 2018
                #fine_slp = interp_finer(reduced_slp2)
                #minslp_lat, minslp_lon, min_slp = min_vals(fine_slp)
                #print [minslp_lat,minslp_lon]
                minslp_lat, minslp_lon, min_slp = find_min_interp(reduced_slp2)
                #print mslp
                #print min_slp
#                print mslp - min_slp
                #exit()
                #print [minslp_lat,minslp_lon]


            if cont_tracking == True and (dd - init_day and hr - init_time > 0)> 3:
                cont_tracking = closed_isobar_check(reduced_slp,minslp_lat, minslp_lon, min_slp)

            max_ws = np.amax(reduced_ws.data)
            min_slp = min_slp/100.
            if cont_tracking == True:
                latitudes.append(minslp_lat); longitudes.append(minslp_lon); minpress.append(min_slp)
                maxws.append(max_ws); valid_days.append(dd); valid_times.append(hr)
            else:
                print "Storm weakened, tracking finished"
                print 'dd = {0:02d}, hr = {1:02d}'.format(dd,hr)
                return latitudes,longitudes, minpress, maxws, valid_days, valid_times
    return latitudes,longitudes, minpress, maxws, valid_days, valid_times


def interp_coords(lat0,lon0,lat1,lon1):
    dlon = lon1 - lon0; dlat = lat1 - lat0
    newlat = lat1 + dlat; newlon = lon1 + dlon
    return newlat,newlon

def max_vals(cube,rad=1e7):
    # Find latitude/longitude coordinates of maximum value of data in cube wihin the defined radius from the centre of the array.
    mid_lat = np.median(cube.coord('latitude').points)
    mid_lon = np.median(cube.coord('longitude').points)
    attempts = 0
    while attempts < 1000:
        index = np.argmax(cube.data)
        indices= np.unravel_index(index,cube.data.shape)
        cenlat = cube.coord('latitude').points[indices[0]] # Find the coordinates for this central position
        cenlon = cube.coord('longitude').points[indices[1]]
        dist = ((cenlat-mid_lat)**2 + (cenlon-mid_lon)**2)**0.5
        if dist < rad:
            maxval = cube.data[indices]
            return cenlat, cenlon, maxval
        else:
            cube.data[indices] = -1e7
        attempts = attempts + 1
    print 'Error: check search radius, ensure grid point is within'
    exit()
    #print maxval


def min_vals(cube,rad=1e7):
    # Find latitude/longitude coordinates of maximum value of data in cube.
    mid_lat = np.median(cube.coord('latitude').points)
    mid_lon = np.median(cube.coord('longitude').points)
    attempts = 0
    while attempts < 1000:
        index = np.argmin(cube.data)
        indices= np.unravel_index(index,cube.data.shape)
        cenlat = cube.coord('latitude').points[indices[0]] # Find the coordinates for this central position
        cenlon = cube.coord('longitude').points[indices[1]]
        dist = ((cenlat-mid_lat)**2 + (cenlon-mid_lon)**2)**0.5
        if dist < rad:
            minval = cube.data[indices]
            return cenlat, cenlon, minval
        else:
            cube.data[indices] = 1e7
        attempts = attempts + 1
    print 'Error: check search radius, ensure grid point is within'
    exit()



def reduce_to_box(minlat,maxlat,minlon,maxlon):
    longitude_constraint1=iris.Constraint(longitude = lambda cell:cell>minlon)
    longitude_constraint2=iris.Constraint(longitude = lambda cell:cell<maxlon)
    latitude_constraint1=iris.Constraint(latitude = lambda cell:cell>minlat)
    latitude_constraint2=iris.Constraint(latitude = lambda cell:cell<maxlat)
    box_constraint=longitude_constraint1&longitude_constraint2&latitude_constraint1&latitude_constraint2
    return box_constraint


def find_init_sc2(md,TT,storm):
    if storm.lower() == 'hagupit':
        yr = 2014
    elif storm.lower() == 'haiyan':
        yr = 2013
    from netCDF4 import Dataset
    from netCDF4 import num2date
    import datetime as dt
    mth = int(str(md)[:2])
    day = int(str(md)[2:])
    fc_date = dt.datetime(yr,mth,day,TT)
    filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm.lower())
    dataset = Dataset(filename)
    lats =  dataset.variables['lat_wmo'][:]
    lons = dataset.variables['lon_wmo'][:]
    times = num2date(dataset.variables['time_wmo'][:],dataset.variables['time_wmo'].units)
    times[:] = [time.replace(microsecond=0) for time in times]

    index = np.where(times == fc_date)
    cenlat = lats[index][0]; cenlon = lons[index][0]
    return cenlat,cenlon

def find_init_sc(md, TT, storm):
    # Find the observed position of the storm on the day and time of the forecast.
    # This will become the first search centre.
    ############ Hagupit ###########################
    # For Hagupit initial date is 30/11/2014, 6Z   #
    ################################################
    ############# Haiyan ############################
    # For Haiyan initial date is 2/11/2013, 6Z      #
    # ***** NEED TO CHECK *****                     #
    #################################################
    if storm =='Hagupit':
        index = 4 * (md - 1201) + (TT / 6) + 3 # i.e. 3 times until 0Z on 1/12, then 4*#days from 1st, then number of 6hr timesteps
    elif storm=='Haiyan':
        index = 4 * (md - 1103) + (TT / 6) + 3
    else:
        print 'Error: Storm not found, check name.'
        exit()
    filename = 'ibtracs/{0}_ibtracs2.nc'.format(storm.lower())
    print 'Fetching track data for Tyhoon {0}'.format(storm)
    dataset = Dataset(filename)
    lats =  dataset.variables['lat_for_mapping'][:]
    lons = dataset.variables['lon_for_mapping'][:]
    cenlat = lats[index]; cenlon = lons[index]
    return cenlat, cenlon


if __name__=='__main__':
main()
