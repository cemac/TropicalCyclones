import iris
import numpy as np
import iris.analysis

def box_constraint(minlat,maxlat,minlon,maxlon):
    # Create constraint to extract data from cube over a certain region
    longitude_constraint1=iris.Constraint(longitude = lambda cell:cell>minlon)
    longitude_constraint2=iris.Constraint(longitude = lambda cell:cell<maxlon)
    latitude_constraint1=iris.Constraint(latitude = lambda cell:cell>minlat)
    latitude_constraint2=iris.Constraint(latitude = lambda cell:cell<maxlat)
    box_constraint=longitude_constraint1&longitude_constraint2&latitude_constraint1&latitude_constraint2
    return box_constraint

def time_constraint(yr,mth,dd,hr):
    time_constraint = iris.Constraint(time=iris.time.PartialDateTime(year=yr,month=mth,day=dd,hour=hr))
    return time_constraint

def get_time(cube):
    # Retrieve time units of cubes
    time_unit = cube.coord('time').units
    validity_time = time_unit.num2date(cube.coord('time').points)[0]
    data_time = time_unit.num2date(cube.coord('forecast_reference_time').points)[0]
    lead_time = cube.coord('forecast_period').points[0]
    return time_unit, validity_time, data_time, lead_time
 
def avg_cubes(*cubes):
    n = float(len(cubes[0]))
    return reduce(iris.analysis.maths.add, cubes[0]) / n    
    
def print_cube(filename):
    cube = iris.load(filename)
    print(cube)
    
def is_hourly_mean(cube):
    # Check to see if data in cube is an hourly mean or not.
    cell_method = iris.coords.CellMethod(method='mean', coords='time',
                                         intervals='1 hour')
    return cell_method in cube.cell_methods


def shift_lon_180to360(cube):
  # Shift longitude range to [-180, 180]
  x_coord = cube.coord(axis="x")
  if np.max(x_coord.points) > 180.0:
    cube = cube.intersection(iris.coords.CoordExtent(x_coord, -180.0, 180.0))
  return cube

def shift_lon_360to180(cube):
  # Shift longitude range to [-180, 180]
  x_coord = cube.coord(axis="x")
  if np.min(x_coord.points) < 0:
    cube = cube.intersection(iris.coords.CoordExtent(x_coord, 0, 360))
  return cube

def average_plevs(cube,ax=0):
    # Averages data along pressure levels taking the full range.
    # Cube should be of the form (pressure,lat,lon) - no time.
    p0 = cube.coord('pressure').points[0]
    p1 = cube.coord('pressure').points[-1]
    dp = p1-p0
    new_data = np.trapz(cube.data,x=cube.coord('pressure').points,axis=ax) / dp
    cube = cube.collapsed('pressure',iris.analysis.MEAN)
    if cube.data.shape != new_data.shape:
        print('Error, check pressure is on axis 0')
        exit()
    cube.data = new_data
    return cube
