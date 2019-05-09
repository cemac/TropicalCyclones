#######################################################################
# TC tracker using UM data. Technique based on that used in the TCVER #
# program by Julian Hemming.                                          #
#                                            John Ashcroft, Sept 2017 #
#######################################################################

import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from netCDF4 import Dataset
import iris
import iris.analysis.calculus
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def tcver(df,lat0,lon0,lat1,lon1):
    plev = 850
    t = 0; dt = 1
    latitudes = []; longitudes = []; min_p = []; max_wspeed = []
    # Constraints for cubes:
    pressure_constraint = iris.Constraint(pressure=plev)
    u_constraint = iris.AttributeConstraint(STASH='m01s15i201')
    v_constraint = iris.AttributeConstraint(STASH='m01s15i202')
    slp_constraint = iris.AttributeConstraint(STASH='m01s16i222')
    sfcu_constraint = iris.AttributeConstraint(STASH='m01s03i225')
    sfcv_constraint = iris.AttributeConstraint(STASH='m01s03i226')

    # Load the cubes
    uwind = iris.load_cube(df,u_constraint).extract(pressure_constraint)
    vwind = iris.load_cube(df,v_constraint).extract(pressure_constraint)
    mslp = iris.load_cube(df,slp_constraint)
    sfc_u = iris.load_cube(df,sfcu_constraint)
    sfc_v = iris.load_cube(df,sfcv_constraint)
    ws = (sfc_u**2 + sfc_v**2)**0.5

    # Find coordintes of domain
    maxlat = np.amax(uwind.coord('latitude').points)
    minlat = np.amin(uwind.coord('latitude').points)
    maxlon = np.amax(uwind.coord('longitude').points)
    minlon = np.amin(uwind.coord('longitude').points)

   # datestimes = datetime.datetime(2014,11,1)


    # Ensure that the time values cover the same points
    print uwind.coord('time').points
    
    exit()
    if len(uwind.coord('time').points) == len(mslp.coord('time').points):
        pass
    elif len(uwind.coord('time').points) == len(mslp.coord('time').points) + 1:
        print 'Changed winds'
        uwind = uwind[1:][:][:]; vwind = vwind[1:][:][:]
    else:
        print 'ERROR: time dimensionts of winds not equal or one greater than time dimensions of slp.'
        exit()

    if len(ws.coord('time').points) == len(mslp.coord('time').points):
        pass
    elif len(ws.coord('time').points) == len(mslp.coord('time').points) + 1:
        ws = ws[1:][:][:];
        print 'Changed windspeeds'
    else:
        print 'ERROR: time dimensions of winds not equal or one greater than time dimensions of slp.'
        exit()
    #print uwind
    # Times are likely to start from T+1, but need to check this.
    for u,v,slp,wspeed in zip(uwind.slices(['latitude','longitude']),vwind.slices(['latitude','longitude']),mslp.slices(['latitude','longitude']),ws.slices(['latitude','longitude'])):
        vrt = iris.analysis.calculus.curl(u,v)[2]
        if t==0:
            box_constraint = reduce_to_box(lat0-3,lat0+3,lon0-3,lon0+3)
            t = t+1
        elif t==1:
            box_constraint = reduce_to_box(lat1-3,lat1+3,lon1-3,lon1+3)
            t = t+1
        else:
            lat2, lon2 = interp_coords(latitudes[-2],longitudes[-2],latitudes[-1],longitudes[-1])
            if lat2 > maxlat or lat2 < minlat or lon2 > maxlon or lon2< minlon:
                print 'WARNING: Storm out of domain'
                return latitudes, longitudes, min_p, max_wspeed
            box_constraint = reduce_to_box(lat2-3,lat2+3,lon2-3,lon2+3)
            t = t+1

        reduced_vrt = vrt.extract(box_constraint)
        v_lat, v_lon, max_vrt = max_vals(reduced_vrt)
        #print 'Max vort position = ({0:.2f},{1:.2f}'.format(v_lat,v_lon)
        box_constraint = reduce_to_box(v_lat-2,v_lat+2,v_lon-2,v_lon+2)
        reduced_slp = slp.extract(box_constraint)
        #print np.amin(reduced_slp.data)
        test_slp_box = reduce_to_box(13,14,127,129)
        test_slp = slp.extract(test_slp_box)
        #print np.amin(test_slp.data)
        c_lat, c_lon, minpress = min_vals(reduced_slp)
        box_constraint = reduce_to_box(c_lat-2,c_lat+2,c_lon-2,c_lon+2)
        reduced_ws = wspeed.extract(box_constraint)
        max_ws = np.amax(reduced_ws.data)
        minpress = minpress / 100.
        print 'T+{0:03d}; lat = {1:.2f}; lon = {2:.2f}; mslp = {3:.0f}; mws = {4:.2f}'.format(dt,c_lat,c_lon,minpress,max_ws)
        dt = dt+6
        latitudes.append(c_lat); longitudes.append(c_lon); min_p.append(minpress); max_wspeed.append(max_ws)
        #print latitudes
        #print min_p
    return latitudes,longitudes, min_p, max_wspeed



def interp_coords(lat0,lon0,lat1,lon1):
    dlon = lon1 - lon0; dlat = lat1 - lat0
    newlat = lat1 + dlat; newlon = lon1 + dlon
    return newlat,newlon

def max_vals(cube):
    # Find latitude/longitude coordinates of maximum value of data in cube.
    index = np.argmax(cube.data)
    indices= np.unravel_index(index,cube.data.shape)
    cenlat = cube.coord('latitude').points[indices[0]] # Find the coordinates for this central position
    cenlon = cube.coord('longitude').points[indices[1]]
    maxval = cube.data[indices]
    #print maxval
    return cenlat, cenlon, maxval

def min_vals(cube):
    # Find latitude/longitude coordinates of maximum value of data in cube.
    index = np.argmin(cube.data)
    indices= np.unravel_index(index,cube.data.shape)
    cenlat = cube.coord('latitude').points[indices[0]] # Find the coordinates for this central position
    cenlon = cube.coord('longitude').points[indices[1]]
    minval = cube.data[indices]
    return cenlat, cenlon, minval


def reduce_to_box(minlat,maxlat,minlon,maxlon):
    longitude_constraint1=iris.Constraint(longitude = lambda cell:cell>minlon)
    longitude_constraint2=iris.Constraint(longitude = lambda cell:cell<maxlon)
    latitude_constraint1=iris.Constraint(latitude = lambda cell:cell>minlat)
    latitude_constraint2=iris.Constraint(latitude = lambda cell:cell<maxlat)
    box_constraint=longitude_constraint1&longitude_constraint2&latitude_constraint1&latitude_constraint2
    return box_constraint

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
    print 'Fetching track data for Tyhoon Hagupit'
    dataset = Dataset(filename)
    lats =  dataset.variables['lat_for_mapping'][:]
    lons = dataset.variables['lon_for_mapping'][:]
    cenlat = lats[index]; cenlon = lons[index]
    return cenlat, cenlon

def main():
############################## Change this #####################################
    monthday = [1203]; times = [12]
    ens_members = np.arange(12)
    models = ['glm','4p4']
    storm = 'Hagupit'
    if storm == 'Hagupit':
        storm_id = 'u-ap087'
    elif storm=='Haiyan':
        storm_id = 'u-ap121'

    for md in monthday:
        for TT in times:
            for mod in models:
                print '***** Finding initial search position. ******'
                lat0, lon0 = find_init_sc(md, TT,storm)
                lat1, lon1 = find_init_sc(md, TT+6,storm)
                print 'Initial search position: ({0:.2f},{1:.2f})'.format(lat0,lon0)
                print '***** Beginning main tcver program. *****'
                for em in ens_members:
                    df = '/nfs/a37/scjea/model_runs/{0}/{1}/data/{2}/{3}_{4:02d}Z/tctracker/tcver_{2}_{3}_{4:02d}Z_em{5:02d}.pp'.format(storm,storm_id,mod,md,TT,em)
                    print 'Ensemble member: {0:02d}'.format(em)
                    lats, lons, minpress, max_ws = tcver(df,lat0,lon0,lat1,lon1)
                while len(lats)<20:
                    lats.append(None); lons.append(None); minpress.append(None); max_ws.append(None)
                np.savez('./track_data/tcver_glm_{0:04d}_{1:02d}Z_en{2:02d}'.format(md,TT,ens_no),lats=lats,lons=lons,minpress=minpress,max_ws=max_ws)
            print '***** End of program *****'

    #md = 1105; TT = 0
    #df = '/nfs/a37/scjea/model_runs/Haiyan/u-ap121/data/4p4/1105_00Z/tctracker/tcver_4p4_1105_00Z_em00.pp'
    #lat0, lon0 = find_init_sc(md, TT)
    #lat1, lon1 = find_init_sc(md, TT+6)
    #lats, lons, minpress, max_ws = tcver(df,lat0,lon0,lat1,lon1)
    #exit()


if __name__=='__main__':
    main()
