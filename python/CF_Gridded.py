# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# Copyright UCAR (c) 2019
# University Corporation for Atmospheric Research(UCAR)
# National Center for Atmospheric Research(NCAR)
# Research Applications Laboratory(RAL)
# P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# 28/01/2019
#
# Name:        module1
# Purpose:
# Author:      $ Kevin Sampson
# Created:     28/01/2019
# Licence:     <your licence>
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

# Written by Kevin Sampson, NCAR
# Modified 2019/01/28
##ALP = True            # If running on Windows machine ALP, must start python VirtualEnv
##if ALP == True:
##    # Activate Virtual Environment so we can use Numpy 1.9.x, which GDAl is compiled against
##    activate_this_file = r'C:\Python27\VirtualEnv_x64\Scripts\activate_this.py'     # Path to Virtual Environment python activate script
##    execfile(activate_this_file, dict(__file__=activate_this_file))

# --- Import Core Modules --- #
import sys
import os
import time
import math

# Import Additional Modules
import numpy
#import netCDF4
import ogr
import osr
# --- End Import Modules --- #

from window_tools import nc_tools


OPT_shift_to_center = True
OPT_add_time_dim = True

# --- Module Configurations --- #
sys.dont_write_bytecode = True                                                  # Do not write compiled (.pyc) files
# --- End Module Configurations --- #

# --- Globals --- #
#BASE_DIR = os.path.join('C:','Users','ksampson','Desktop','WINDOW_data','CF_Converter')
#if not os.path.exists(BASE_DIR):
#    BASE_DIR = os.getcwd()
#in_nc = os.path.join(BASE_DIR, 'ncar_les_geo_em.d01.nc')
#projdir = BASE_DIR
#out_nc = os.path.join(projdir, os.path.basename(in_nc).replace('.nc', '_spatial.nc'))
geogridVariable = 'HGT_M'                                                       # 2D Variable in the GEOGRID file to use in defining the grid
processing_notes_SM = 'Created for testing WINDOW AIR converter'                # Notes section appended to the output netCDF file as a global attribute
XYmap = {'x': 'west_east', 'y': 'south_north'}                                  # Dictionary to define mapping between (x,y) and input file dimension names

# Latitude and Longitude 2D array information
addLatLon = True                                                                # Add Latitude and Longitude variables to output file?
addData = True                                                                  # Switch to add data variable defined in 'geogridVariable' to output
latVar = 'XLAT_M'                                                               # Latitude variable in the input file (must be on grid cell center)
lonVar = 'XLONG_M'                                                              # Longitude variable in the input file (must be on grid cell center)

# Initiate dictionaries of GEOGRID projections and parameters
#   See http://www.mmm.ucar.edu/wrf/users/docs/user_guide_V3/users_guide_chap3.htm#_Description_of_the_1
projdict = {1: 'Lambert Conformal Conic',
            2: 'Polar Stereographic',
            3: 'Mercator',
            6: 'Cylindrical Equidistant'}
CF_projdict = {1: "lambert_conformal_conic",
               2: "polar_stereographic",
               3: "mercator",
               6: "latitude_longitude",
               0: "crs"}

# Unify all coordinate system variables to have the same name ("crs"). Ths makes it easier for WRF-Hydro output routines to identify the variable and transpose it to output files
crsVarname = True                                                               # Switch to make all coordinate system variables = "crs" instead of related to the coordinate system name
crsVar = CF_projdict[0]                                                         # Expose this as a global for other functions in other scripts to use
CFConv = 'CF-1.5'                                                               # CF-Conventions version to place in the 'Conventions' attribute of RouteLink files. Maybe 1.0 is enough?
PpVersion = 'v1 (01/2019)'                                                      # WINDOW AIR version to add to output metadata
outNCType = 'NETCDF4_CLASSIC'                                                   # Define the output netCDF version for RouteLink.nc and LAKEPARM.nc

# Global attributes for altering the sphere radius used in computations. Do not alter sphere_radius for standard WRF simulations
sphere_radius = 6370000.0                                                       # Radius of sphere to use (WRF Default = 6370000.0m)
# --- End Globals --- #

# --- Functions --- #
def Read_GEOGRID_for_SRS(in_nc):
    """
    The input NetCDF file (Ideally WRF GEOGRID) gets georeferenced and projection
    infromation is created.

    10/6/2017: Proj4 string generation was added, with definitions adapted from
    https://github.com/NCAR/wrf-python/blob/develop/src/wrf/projection.py
    """

    tic = time.time()

    # First step: Import and georeference NetCDF file
    print('  Step 1: NetCDF Conversion initiated...')
    print('    Input netCDF GEOGRID file: %s' %in_nc)

    # Read input WPS GEOGRID file
    # Loop through global variables in NetCDF file to gather projection information
    rootgrp = nc_tools.open_nc_Dataset(in_nc)                                   # Establish an object for reading the input NetCDF file
    globalAtts = rootgrp.__dict__                                               # Read all global attributes into a dictionary
    map_pro = globalAtts['MAP_PROJ']                                            # Find out which projection this GEOGRID file is in
    print('    Map Projection: %s' %projdict[map_pro])

    # Collect grid corner XY and DX DY for creating ascii raster later
    if 'corner_lats' in globalAtts:
        corner_lat = globalAtts['corner_lats'][13].astype(numpy.float64)        # Note: The values returned are corner points of the mass grid. 13 = Upper left of the Unstaggered grid
    if 'corner_lons' in globalAtts:
        corner_lon = globalAtts['corner_lons'][13].astype(numpy.float64)        # Note: The values returned are corner points of the mass grid. 13 = Upper left of the Unstaggered grid
    if 'DX' in globalAtts:
        DX = globalAtts['DX'].astype(numpy.float32)
    if 'DY' in globalAtts:
        DY = globalAtts['DY'].astype(numpy.float32)

    # Collect necessary information to put together the projection file
    if 'TRUELAT1' in globalAtts:
        standard_parallel_1 = globalAtts['TRUELAT1'].astype(numpy.float64)
    if 'TRUELAT2' in globalAtts:
        standard_parallel_2 = globalAtts['TRUELAT2'].astype(numpy.float64)
    if 'STAND_LON' in globalAtts:
        central_meridian = globalAtts['STAND_LON'].astype(numpy.float64)
    if 'POLE_LAT' in globalAtts:
        pole_latitude = globalAtts['POLE_LAT'].astype(numpy.float64)
    if 'POLE_LON' in globalAtts:
        pole_longitude = globalAtts['POLE_LON'].astype(numpy.float64)
    if 'MOAD_CEN_LAT' in globalAtts:
        print('    Using MOAD_CEN_LAT for latitude of origin.')
        latitude_of_origin = globalAtts['MOAD_CEN_LAT'].astype(numpy.float64)         # Added 2/26/2017 by KMS
    elif 'CEN_LAT' in globalAtts:
        print('    Using CEN_LAT for latitude of origin.')
        latitude_of_origin = globalAtts['CEN_LAT'].astype(numpy.float64)
    del globalAtts
    rootgrp.close()

    # Initiate OSR spatial reference object - See http://gdal.org/java/org/gdal/osr/SpatialReference.html
    proj1 = osr.SpatialReference()

    # Use projection information from global attributes to populate OSR spatial reference object
    # See this website for more information on defining coordinate systems: http://gdal.org/java/org/gdal/osr/SpatialReference.html
    print('    Map Projection: %s' %projdict[int(map_pro)])
    if map_pro == 1:
        # Lambert Conformal Conic
        if 'standard_parallel_2' in locals():
            proj1.SetLCC(standard_parallel_1, standard_parallel_2, latitude_of_origin, central_meridian, 0, 0)
            #proj1.SetLCC(double stdp1, double stdp2, double clat, double clong, double fe, double fn)        # fe = False Easting, fn = False Northing
        else:
            proj1.SetLCC1SP(latitude_of_origin, central_meridian, 1, 0, 0)       # Scale = 1???
            #proj1.SetLCC1SP(double clat, double clong, double scale, double fe, double fn)       # 1 standard parallell
    elif map_pro == 2:
        # Polar Stereographic
        phi1 = standard_parallel_1
        ### Back out the central_scale_factor (minimum scale factor?) using formula below using Snyder 1987 p.157 (USGS Paper 1395)
        ##phi = math.copysign(float(pole_latitude), float(latitude_of_origin))    # Get the sign right for the pole using sign of CEN_LAT (latitude_of_origin)
        ##central_scale_factor = (1 + (math.sin(math.radians(phi1))*math.sin(math.radians(phi))) + (math.cos(math.radians(float(phi1)))*math.cos(math.radians(phi))))/2
        # Method where central scale factor is k0, Derivation from C. Rollins 2011, equation 1: http://earth-info.nga.mil/GandG/coordsys/polar_stereographic/Polar_Stereo_phi1_from_k0_memo.pdf
        # Using Rollins 2011 to perform central scale factor calculations. For a sphere, the equation collapses to be much  more compact (e=0, k90=1)
        central_scale_factor = (1 + math.sin(math.radians(abs(phi1))))/2        # Equation for k0, assumes k90 = 1, e=0. This is a sphere, so no flattening
        print('        Central Scale Factor: %s' %central_scale_factor)
        #proj1.SetPS(latitude_of_origin, central_meridian, central_scale_factor, 0, 0)    # example: proj1.SetPS(90, -1.5, 1, 0, 0)
        proj1.SetPS(pole_latitude, central_meridian, central_scale_factor, 0, 0)    # Adjusted 8/7/2017 based on changes made 4/4/2017 as a result of Monaghan's polar sterographic domain. Example: proj1.SetPS(90, -1.5, 1, 0, 0)
        #proj1.SetPS(double clat, double clong, double scale, double fe, double fn)
    elif map_pro == 3:
        # Mercator Projection
        proj1.SetMercator(latitude_of_origin, central_meridian, 1, 0, 0)     # Scale = 1???
        #proj1.SetMercator(double clat, double clong, double scale, double fe, double fn)
    elif map_pro == 6:
        # Cylindrical Equidistant (or Rotated Pole)
        if pole_latitude != float(90) or pole_longitude != float(0):
            # if pole_latitude, pole_longitude, or stand_lon are changed from thier default values, the pole is 'rotated'.
            print('[PROBLEM!] Cylindrical Equidistant projection with a rotated pole is not currently supported.')
            raise SystemExit
        else:
            proj1.SetEquirectangular(latitude_of_origin, central_meridian, 0, 0)
            #proj1.SetEquirectangular(double clat, double clong, double fe, double fn)
            #proj1.SetEquirectangular2(double clat, double clong, double pseudostdparallellat, double fe, double fn)

    # Set Geographic Coordinate system (datum) for projection
    proj1.SetGeogCS('WRF_Sphere', 'Sphere', '', sphere_radius, 0.0)      # Could try 104128 (EMEP Sphere) well-known?
    #proj1.SetGeogCS(String pszGeogName, String pszDatumName, String pszSpheroidName, double dfSemiMajor, double dfInvFlattening)

    # Set the origin for the output raster (in GDAL, usuall upper left corner) using projected corner coordinates
    wgs84_proj = osr.SpatialReference()
    wgs84_proj.ImportFromEPSG(4326)
    transform = osr.CoordinateTransformation(wgs84_proj, proj1)
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint_2D(corner_lon, corner_lat)
    point.Transform(transform)
    x00 = point.GetX(0)
    y00 = point.GetY(0)
    print('  Created projection definition from input NetCDF GEOGRID file %s in %.2fs.' %(in_nc, time.time()-tic))
    return proj1, DX, DY, x00, y00, map_pro

def add_CRS_var(rootgrp, sr, map_pro, CoordSysVarName, grid_mapping, PE_stringEsri, PE_stringGDAL,
                GeoTransformStr=None, coordinateAxesName='y x'):
    """
    This function was added to generalize the creating of a coordinate reference
    system variable based on a spatial reference object and other grid and projection
    information.

    To find OGR projection parameters, reference https://www.gdal.org/ogr__srs__api_8h.html
    """
    tic1 = time.time()

    # Scalar projection variable - http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
    proj_var = rootgrp.createVariable(CoordSysVarName, 'S1')                    # (Scalar Char variable)
    proj_var.transform_name = grid_mapping                                      # grid_mapping. grid_mapping_name is an alias for this
    proj_var.grid_mapping_name = grid_mapping                                   # for CF compatibility
    proj_var.esri_pe_string = PE_stringEsri                                     # For ArcGIS. Not required if esri_pe_string exists in the 2D variable attributes
    proj_var.spatial_ref = PE_stringGDAL                                        # For GDAl
    proj_var.long_name = "CRS definition"                                       # Added 10/13/2017 by KMS to match GDAL format
    proj_var.longitude_of_prime_meridian = 0.0                                  # Added 10/13/2017 by KMS to match GDAL format
    if GeoTransformStr is not None:
        proj_var.GeoTransform = GeoTransformStr                                 # For GDAl - GeoTransform array

    # Projection specific parameters - http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
    if map_pro == 1:
        # Lambert Conformal Conic

        # Required transform variables
        proj_var._CoordinateAxes = coordinateAxesName                           # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        proj_var._CoordinateTransformType = "Projection"
        proj_var.standard_parallel = sr.GetProjParm("standard_parallel_1"), sr.GetProjParm("standard_parallel_2")     # Double
        proj_var.longitude_of_central_meridian = float(sr.GetProjParm("central_meridian")) # Double. Necessary in combination with longitude_of_prime_meridian?
        proj_var.latitude_of_projection_origin = float(sr.GetProjParm("latitude_of_origin"))         # Double

        # Optional tansform variable attributes
        proj_var.false_easting = float(sr.GetProjParm("false_easting"))         # Double  Always in the units of the x and y projection coordinates
        proj_var.false_northing = float(sr.GetProjParm("false_northing"))       # Double  Always in the units of the x and y projection coordinates
        proj_var.earth_radius = sphere_radius                                   # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km.
        proj_var.semi_major_axis = sphere_radius                                # Added 10/13/2017 by KMS to match GDAL format
        proj_var.inverse_flattening = float(0)                                  # Added 10/13/2017 by KMS to match GDAL format: Double - optional Lambert Conformal Conic parameter

    elif map_pro == 2:
        # Polar Stereographic

        # Required transform variables
        proj_var._CoordinateAxes = coordinateAxesName                           # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        proj_var._CoordinateTransformType = "Projection"
        proj_var.longitude_of_projection_origin = float(sr.GetProjParm("longitude_of_origin"))   # Double - proj_var.straight_vertical_longitude_from_pole = ''
        proj_var.latitude_of_projection_origin = float(sr.GetProjParm("latitude_of_origin"))     # Double
        proj_var.scale_factor_at_projection_origin = float(sr.GetProjParm("scale_factor"))      # Double

        # Optional tansform variable attributes
        proj_var.false_easting = float(sr.GetProjParm("false_easting"))                         # Double  Always in the units of the x and y projection coordinates
        proj_var.false_northing = float(sr.GetProjParm("false_northing"))                       # Double  Always in the units of the x and y projection coordinates
        proj_var.earth_radius = sphere_radius                                   # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km.
        proj_var.semi_major_axis = sphere_radius                                # Added 10/13/2017 by KMS to match GDAL format
        proj_var.inverse_flattening = float(0)                                  # Added 10/13/2017 by KMS to match GDAL format: Double - optional Lambert Conformal Conic parameter

    elif map_pro == 3:
        # Mercator

        # Required transform variables
        proj_var._CoordinateAxes = coordinateAxesName                           # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        proj_var._CoordinateTransformType = "Projection"
        proj_var.longitude_of_projection_origin = float(sr.GetProjParm("longitude_of_origin"))   # Double
        proj_var.latitude_of_projection_origin = float(sr.GetProjParm("latitude_of_origin"))     # Double
        proj_var.standard_parallel = float(sr.GetProjParm("standard_parallel_1"))                # Double
        proj_var.earth_radius = sphere_radius                                   # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km.
        proj_var.semi_major_axis = sphere_radius                                # Added 10/13/2017 by KMS to match GDAL format
        proj_var.inverse_flattening = float(0)                                  # Added 10/13/2017 by KMS to match GDAL format: Double - optional Lambert Conformal Conic parameter

    elif map_pro == 6:
        # Cylindrical Equidistant or rotated pole

        #http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#appendix-grid-mappings
        # Required transform variables
        #proj_var.grid_mapping_name = "latitude_longitude"                      # or "rotated_latitude_longitude"

        print('        Cylindrical Equidistant projection not supported.')
        raise SystemExit

    # Added 10/13/2017 by KMS to accomodate alternate datums
    elif map_pro == 0:
        proj_var._CoordinateAxes = 'lat lon'
        proj_var.semi_major_axis = sr.GetSemiMajor()
        proj_var.semi_minor_axis =  sr.GetSemiMinor()
        if sr.flattening != 0:
            proj_var.inverse_flattening = float(float(1)/sr.GetInvFlattening()) # This avoids a division by 0 error
        else:
            proj_var.inverse_flattening = float(0)
        pass

    if OPT_add_time_dim:
        postfix ='_t'
        proj_t_var = nc_tools.copy_variable(rootgrp, proj_var, var_name=CoordSysVarName + postfix)
        #proj_t_var.transform_name = grid_mapping + postfix                          # grid_mapping. grid_mapping_name is an alias for this
        proj_t_var.transform_name = grid_mapping                                     # grid_mapping. grid_mapping_name is an alias for this
        proj_t_var.grid_mapping_name = grid_mapping + postfix                        # for CF compatibility
        proj_t_var.long_name = "CRS definition with time"                            # Added 10/13/2017 by KMS to match GDAL format
        #proj_t_var._CoordinateAxes = 't ' + coordinateAxesName                        # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        proj_t_var._CoordinateAxes = 't ' + coordinateAxesName                        # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        proj_t_var._CoordinateAxes = 'XTIME ' + coordinateAxesName                    # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems

    # Global attributes related to CF-netCDF
    rootgrp.Conventions = CFConv
    return rootgrp

def getXY(GeoTransformStr, nx, ny, flipY=True):
    """
    This function will use the affine transformation (GeoTransform) to produce an
    array of X and Y 1D arrays. Note that the GDAL affine transformation provides
    the grid cell coordinates from the upper left corner. This is typical in GIS
    applications. However, WRF uses a south_north ordering, where the arrays are
    written from the bottom to the top. Thus, in order to flip the y array, select
    flipY = True (default).
    """
    tic1 = time.time()
    gt = [float(item) for item in GeoTransformStr.split(' ')]

    # Build i,j arrays
    #j = numpy.arange(dataArr.shape[0])
    #i = numpy.arange(dataArr.shape[1])
    j = numpy.arange(nx)
    i = numpy.arange(ny)
    if OPT_shift_to_center:
        j = j + float(0.5)                             # Add 0.5 to estimate coordinate of grid cell centers
        i = i + float(0.5)                             # Add 0.5 to estimate coordinate of grid cell centers

    # col, row to x, y   From https://www.perrygeo.com/python-affine-transforms.html
    x = (i * gt[1]) + gt[0]
    y = (j * gt[5]) + gt[3]

    if flipY:
        y = y[::-1]
    print('    1D Coordinate variables calculated in {0: 8.2f} seconds.'.format(time.time()-tic1))
    #print('     DEBUG ====  getXY: x: ', x[0] , ' to  ', x[-1], ' y: ', y[0],' to ', y[-1], ' length=', len(y))
    #print('     DEBUG ====  getXY: from gt[1]: ', gt[1] , ' gt[0]: ', gt[0], ' gt[5]: ', gt[5],' gt[3]: ', gt[3])
    return x, y


def create_CF_NetCDF(rootgrp, dataArr, sr, map_pro, projdir, DX, DY, GeoTransformStr,
                     Xsize, Ysize, addLatLon=False, addVars=[],
                     WKT_Esri="", WKT_GDAL="", addData=False, GeoTransform_alt=None,
                     Xsize_alt=-1, Ysize_alt=-1):
    """This function will create the netCDF file with CF conventions for the grid
    description. The output NetCDF will have the XMAP/YMAP created for the x and
    y variables and the LATITUDE and LONGITUDE variables populated with the selected
    global variables latVar and lonVar from the GEOGRID file."""

    tic1 = time.time()
    print('  Creating CF-netCDF File.')

    # Find name for the grid mapping
    if CF_projdict.get(map_pro) is not None:
        grid_mapping = CF_projdict[map_pro]
        print('    Map Projection of input raster : %s' %grid_mapping)
    else:
        grid_mapping = 'crs'                                                    # Added 10/13/2017 by KMS to generalize the coordinate system variable names
        print('    Map Projection of input raster (not a WRF projection): %s' %grid_mapping)

    # Create Dimensions
    reg_dx = float(DX)
    reg_dy = float(DY)
    sub_dx = reg_dx
    sub_dy = reg_dy
    dim_y = rootgrp.createDimension('y', Ysize)
    dim_x = rootgrp.createDimension('x', Xsize)
    has_subgrid = (0 < Ysize_alt) and (0 < Xsize_alt)
    if has_subgrid:
        my_subgrid_factor = int(round(Xsize_alt / Xsize)) if Xsize_alt > Xsize else int(round(Xsize / Xsize_alt)) 
        dim_y_sub = rootgrp.createDimension('y_alt', Ysize_alt)
        dim_x_sub = rootgrp.createDimension('x_alt', Xsize_alt)
        if Xsize > Xsize_alt:
            reg_dx = reg_dx / my_subgrid_factor
        else:
            sub_dx = reg_dx / my_subgrid_factor
        if Ysize > Ysize_alt:
            reg_dy = reg_dy / my_subgrid_factor
        else:
            sub_dy = reg_dy / my_subgrid_factor
    print('    Dimensions created after {0: 8.2f} seconds.'.format(time.time()-tic1))

    append_vars = []

    # Create coordinate variables
    #var_y = rootgrp.createVariable('y', 'f8', 'y')                              # (64-bit floating point)
    #var_x = rootgrp.createVariable('x', 'f8', 'x')                              # (64-bit floating point)
    var_y = rootgrp.createVariable('y', 'f8', dim_y.name)                        # (64-bit floating point)
    var_x = rootgrp.createVariable('x', 'f8', dim_x.name)                        # (64-bit floating point)
    append_vars.append('x')
    append_vars.append('y')
    if has_subgrid:
        var_y_subgrid = rootgrp.createVariable('y_alt', 'f8', dim_y_sub.name)        # (64-bit floating point)
        var_x_subgrid = rootgrp.createVariable('x_alt', 'f8', dim_x_sub.name)        # (64-bit floating point)
        append_vars.append('x_alt')
        append_vars.append('y_alt')

    # Must handle difference between ProjectionCoordinateSystem and LatLonCoordinateSystem
    if proj.IsGeocentric():
        if crsVarname:
            CoordSysVarName = crsVar
        else:
            CoordSysVarName = "LatLonCoordinateSystem"

        # Set variable attributes
        #var_y.standard_name = ''
        #var_x.standard_name = ''
        var_y.long_name = "latitude coordinate"
        var_x.long_name = "longitude coordinate"
        var_y.units = "degrees_north"
        var_x.units = "degrees_east"
        var_y._CoordinateAxisType = "Lat"
        var_x._CoordinateAxisType = "Lon"
        if var_y_subgrid is not None:
            var_y_subgrid.long_name = "latitude sub coordinate"
            var_y_subgrid.units = var_y.units
            var_y_subgrid._CoordinateAxisType = var_y._CoordinateAxisType
        if var_x_subgrid is not None:
            var_x_subgrid.long_name = "longitude sub coordinate"
            var_x_subgrid.units = var_x.units
            var_x_subgrid._CoordinateAxisType = var_x._CoordinateAxisType

    elif proj.IsProjected():
        if crsVarname:
            CoordSysVarName = crsVar
        else:
            CoordSysVarName = "ProjectionCoordinateSystem"
        #proj_units = sr.linearUnitName.lower()                                  # sr.projectionName wouldn't work for a GEOGCS
        proj_units = 'm'                                                        # Change made 11/3/2016 by request of NWC

        # Set variable attributes
        var_y.standard_name = 'projection_y_coordinate'
        var_x.standard_name = 'projection_x_coordinate'
        var_y.long_name = 'y coordinate of projection'
        var_x.long_name = 'x coordinate of projection'
        var_y.units = proj_units                                                # was 'meter', now 'm'
        var_x.units = proj_units                                                # was 'meter', now 'm'
        var_y._CoordinateAxisType = "GeoY"                                      # Use GeoX and GeoY for projected coordinate systems only
        var_x._CoordinateAxisType = "GeoX"                                      # Use GeoX and GeoY for projected coordinate systems only
        var_y.resolution = reg_dy                                               # Added 11/3/2016 by request of NWC
        var_x.resolution = reg_dx                                               # Added 11/3/2016 by request of NWC

        if var_y_subgrid is not None:
            var_y_subgrid.standard_name = 'projection_y_alt_coordinate'
            var_y_subgrid.long_name = 'y_alt coordinate of projection'
            var_y_subgrid.units = var_y.units
            var_y_subgrid._CoordinateAxisType = var_y._CoordinateAxisType
            var_y_subgrid.resolution = sub_dy
        if var_x_subgrid is not None:
            var_x_subgrid.standard_name = 'projection_x_alt_coordinate'
            var_x_subgrid.long_name = 'x_alt coordinate of projection'
            var_x_subgrid.units = var_x.units
            var_x_subgrid._CoordinateAxisType = var_x._CoordinateAxisType
            var_x_subgrid.resolution = sub_dx

        # Build coordinate reference system variable
        rootgrp = add_CRS_var(rootgrp, sr, map_pro, CoordSysVarName, grid_mapping, WKT_Esri, WKT_GDAL, GeoTransformStr)
        append_vars.append(CoordSysVarName)
        if rootgrp.variables.get(CoordSysVarName + "_t", None) is not None:
            append_vars.append(CoordSysVarName + "_t")
        if GeoTransform_alt is not None:
            sub_CoordSysVarName = CoordSysVarName + "_alt"
            add_CRS_var(rootgrp, sr, map_pro, sub_CoordSysVarName, grid_mapping, WKT_Esri, WKT_GDAL, GeoTransform_alt, 'y_alt x_alt')
            append_vars.append(sub_CoordSysVarName)
            rootgrp.setncattr('WINDOW_grid_mapping_alt', sub_CoordSysVarName)
            rootgrp.setncattr('WINDOW_CoordinateSystems_alt', sub_CoordSysVarName)
            if rootgrp.variables.get(sub_CoordSysVarName + "_t", None) is not None:
                append_vars.append(sub_CoordSysVarName + "_t")

    if 0 < len(WKT_Esri):
        rootgrp.setncattr('WINDOW_esri_pe_string', WKT_Esri)
    rootgrp.setncattr('WINDOW_grid_mapping', CoordSysVarName)
    rootgrp.setncattr('WINDOW_CoordinateSystems', CoordSysVarName)

    # For prefilling additional variables and attributes on the same 2D grid, given as a list [[<varname>, <vardtype>, <long_name>],]
    for varinfo in addVars:
        ncvar = rootgrp.createVariable(varinfo[0], varinfo[1], ('y', 'x'))
        ncvar.esri_pe_string = WKT_Esri
        ncvar.grid_mapping = CoordSysVarName
        #ncvar.long_name = varinfo[2]
        #ncvar.units = varinfo[3]

    print('    Coordinate variables and variable attributes set after {0: 8.2f} seconds.'.format(time.time()-tic1))

    dim_names = ('y', 'x')
    if Xsize > Xsize_alt:
        dim_names = ('y_alt', 'x_alt')
    if addLatLon == True:

        print('    Proceeding to add LATITUDE and LONGITUDE variables after {0: 8.2f} seconds.'.format(time.time()-tic1))

        # Populate this file with 2D latitude and longitude variables
        # Latitude and Longitude variables (WRF)
        lat_WRF = rootgrp.createVariable('LATITUDE', 'f4', dim_names)           # (32-bit floating point)
        lon_WRF = rootgrp.createVariable('LONGITUDE', 'f4', dim_names)          # (32-bit floating point)
        lat_WRF.long_name = 'latitude coordinate'                               # 'LATITUDE on the WRF Sphere'
        lon_WRF.long_name = 'longitude coordinate'                              # 'LONGITUDE on the WRF Sphere'
        lat_WRF.units = "degrees_north"
        lon_WRF.units = "degrees_east"
        lat_WRF._CoordinateAxisType = "Lat"
        lon_WRF._CoordinateAxisType = "Lon"
        lat_WRF.grid_mapping = CoordSysVarName                                  # This attribute appears to be important to Esri
        lon_WRF.grid_mapping = CoordSysVarName                                  # This attribute appears to be important to Esri

        '''Adding the Esri PE String in addition to the CF grid mapping attributes
        is very useful. Esri will prefer the PE string over other CF attributes,
        allowing a spherical datum to be defined. Esri can interpret the coordinate
        system variable alone, but will assume the datum is WGS84. This cannot be
        changed except when using an Esri PE String.'''

        lat_WRF.esri_pe_string = WKT_Esri
        lon_WRF.esri_pe_string = WKT_Esri

        # Missing value attribute not needed yet
        #missing_val = numpy.finfo(numpy.float32).min                            # Define missing data variable based on numpy
        #lat_WRF.missing_value = missing_val                                     # Float sys.float_info.min?
        #lon_WRF.missing_value = missing_val                                     # Float sys.float_info.min?

        ##    # Create a new coordinate system variable
        ##    LatLonCoordSysVarName = "LatLonCoordinateSystem"
        ##    latlon_var = rootgrp.createVariable(LatLonCoordSysVarName, 'S1')            # (Scalar Char variable)
        ##    latlon_var._CoordinateAxes = 'LATITUDE LONGITUDE'                           # Coordinate systems variables always have a _CoordinateAxes attribute

        # Data variables need _CoodinateSystems attribute
        lat_WRF._CoordinateAxisType = "Lat"
        lon_WRF._CoordinateAxisType = "Lon"
        lat_WRF._CoordinateSystems = CoordSysVarName
        lon_WRF._CoordinateSystems = CoordSysVarName
        ##    lat_WRF._CoordinateSystems = "%s %s" %(CoordSysVarName, LatLonCoordSysVarName)        # For specifying more than one coordinate system
        ##    lon_WRF._CoordinateSystems = "%s %s" %(CoordSysVarName, LatLonCoordSysVarName)        # For specifying more than one coordinate system

    if addData:
        data_Var = rootgrp.createVariable(geogridVariable, dataArr.dtype, dim_names)          # (32-bit floating point)
        data_Var.grid_mapping = CoordSysVarName                                 # This attribute appears to be important to Esri
        data_Var.esri_pe_string = WKT_Esri
        data_Var._CoordinateSystems = CoordSysVarName
        #data_Var.units = ''
        #data_Var.long_name = ''
        data_Var[:] = dataArr
        #append_vars.append(geogridVariable)


    print('    netCDF global attributes set after {0: 8.2f} seconds.'.format(time.time()-tic1))
    return rootgrp, append_vars

try:
    #from osgeo import ogr, osr, gdal
    from osgeo import gdal
except:
    sys.exit('ERROR: cannot find GDAL/OGR modules')
    
# example GDAL error handler function
def gdal_error_handler(err_class, err_num, err_msg):
    errtype = {
            gdal.CE_None:'None',
            gdal.CE_Debug:'Debug',
            gdal.CE_Warning:'Warning',
            gdal.CE_Failure:'Failure',
            gdal.CE_Fatal:'Fatal'
    }
    err_msg = err_msg.replace('\n',' ')
    err_class = errtype.get(err_class, 'None')
    print('Error Number: %s' % (err_num))
    print('Error Type: %s' % (err_class))
    print('Error Message: %s' % (err_msg))
    
def process_arguments(argv):
    arg_index = 0
    arg_count = len(argv)
    in_nc = argv[arg_index]
    if not os.path.exists(in_nc):
        print(" ==ERROR == The input geogrid file [{i}] does not exist!!!".format(i=in_nc))
        print("            Quit...")
        sys.exit(-2);
    
    out_nc = None
    sub_grid_factor = 1
    subgrid_first = False
    
    arg_index += 1
    projdir = os.path.dirname(in_nc)
    if 2 < arg_count:
        if not argv[arg_index].isdigit() and not argv[arg_index].startswith('subgrid_'):
            out_nc = argv[arg_index]
            arg_index += 1
        if argv[arg_index].isdigit():
            sub_grid_factor = int(argv[arg_index])
            arg_index += 1
        elif argv[arg_index].startswith('subgrid_factor'):
            sub_grid_factor = int(argv[arg_index].split('=')[1])
            arg_index += 1
        if (arg_index < arg_count) and 'subgrid_first' == argv[arg_index]:
            subgrid_first = True
    if out_nc is None:
        out_nc = os.path.join(projdir, os.path.basename(in_nc).replace('.nc', '_spatial.nc'))
    return in_nc, out_nc, sub_grid_factor, subgrid_first

def show_help(app_name):
    print(" Usage {p} input_geogrid_file [output_geogrid_file]".format(p=app_name))
    print("     input_geogrid_file: input geogrid file name, required")
    print("    output_geogrid_file: output geogrid file name, optional, default: <input>_spatial.nc")
    print("   <subgrid_factor=DD>: integer to make sub_grid dimension, optional")
    print("   <subgrid_first>: string 'subgrid_first' to make sub_grid dimension, optional")
    
if __name__ == '__main__':
    debug = False
    #debug = not debug
    tic = time.time()

    # Gather all necessary parameters
    if 1 == len(sys.argv):
        print(" Usage {p} input_geogrid_file [output_geogrid_file]".format(p=sys.argv[0]))
        print("     input_geogrid_file: input geogrid file name, required")
        print("    output_geogrid_file: output geogrid file name, optional, default: <input>_spatial.nc")
        print("   <sub_grid_factor=DD>: integer to make sub_grid dimension, optional")
        sys.exit(-1);
    
    in_nc, out_nc, sub_grid_factor, subgrid_first = process_arguments(sys.argv[1:])
    if debug:
        print(' DEBUG: sub_grid_factor: {f} subgrid_first: {s}'.format(f=sub_grid_factor, s=subgrid_first))
    projdir = os.path.dirname(in_nc)

    # install error handler
    gdal.PushErrorHandler(gdal_error_handler)
    
    # Get projection information as an OSR SpatialReference object from Geogrid File
    proj, DX, DY, x00, y00, map_pro = Read_GEOGRID_for_SRS(in_nc)

    # Build coordinate system representations
    proj4 = proj.ExportToProj4()                                                # Proj.4 string coordinate system representation
    GeoTransform = '%s %s %s %s %s %s' %(x00, DX, 0, y00, 0, -DY)               # Build an affine transformation (useful info to add to metadata)
    PE_stringGDAL = proj.ExportToWkt()                                          # This is GDAL's representation of a WTK string
    proj.MorphToESRI()                                                          # Alter the projection to Esri's representation of a coordinate system
    PE_stringEsri = proj.ExportToWkt()                                          # This is Esri's representation of a WKT string
    print('  Proj4: {0}'.format(proj4))                                         # Print Proj.4 string to screen
    print('  GeoTransform: {0}'.format(GeoTransform))                           # Print affine transformation to screen.
    if sub_grid_factor > 1:
        GeoTransform_subgrid = '%s %s %s %s %s %s' %(x00, DX/sub_grid_factor, 0, y00, 0, -DY/sub_grid_factor)               # Build an affine transformation (useful info to add to metadata)

    # Gather variable from the input file
    rootgrp_in = nc_tools.open_nc_Dataset(in_nc)                                # Open read object on input netCDF file
    ncvar = rootgrp_in.variables[geogridVariable]                               # netCDF Variable object
    varDims = ncvar.dimensions                                                  # Read variable dimensions for the selected GEOGRID variable
    Ysize = ncvar.shape[ncvar.dimensions.index(XYmap['y'])]                     # Get the dimension size based on the dimension name
    Xsize = ncvar.shape[ncvar.dimensions.index(XYmap['x'])]                     # Get the dimension size based on the dimension name
    dataArr = ncvar[0]                                                          # Select first time slice. Assumes time is first dimension.
    del ncvar

    Xsize_reg = Xsize
    Ysize_reg = Ysize
    has_subgrid = sub_grid_factor > 1
    GeoTransform_reg = GeoTransform
    GeoTransform_alt = None
    if has_subgrid:
        GeoTransform_alt = GeoTransform_subgrid
        Xsize_subgrid = (Xsize+1)*sub_grid_factor
        Ysize_subgrid = (Ysize+1)*sub_grid_factor
        if subgrid_first:
            Xsize_alt = Xsize
            Ysize_alt = Ysize
            Xsize_reg = Xsize_subgrid
            Ysize_reg = Ysize_subgrid
            GeoTransform_alt = GeoTransform
            GeoTransform_reg = GeoTransform_subgrid
        else:
            Xsize_alt = Xsize_subgrid
            Ysize_alt = Ysize_subgrid
    # Gather 1D coordinates using data array size and geotransform
    xarr, yarr = getXY(GeoTransform_reg, Xsize_reg, Ysize_reg, flipY=True)                       # Use affine transformation to calculate cell coordinates

    # Create spatial metadata file for GEOGRID/LDASOUT grids
    rootgrp = nc_tools.create_nc_Dataset(out_nc, nc_format=outNCType)                    # Open write object (create) on output netCDF file
    rootgrp, append_vars = create_CF_NetCDF(rootgrp, dataArr, proj, map_pro, projdir, DX, DY,
        GeoTransform_reg, Xsize_reg, Ysize_reg, addLatLon=addLatLon, WKT_Esri=PE_stringEsri,
        WKT_GDAL=PE_stringGDAL, addData=addData, GeoTransform_alt=GeoTransform_alt,
        Xsize_alt=Xsize_alt, Ysize_alt=Ysize_alt)

    if addLatLon:
        rootgrp.variables['LATITUDE'][:] = rootgrp_in.variables[latVar][0]
        rootgrp.variables['LONGITUDE'][:] = rootgrp_in.variables[lonVar][0]
        append_vars.append('LATITUDE')
        append_vars.append('LONGITUDE')
    rootgrp_in.close()

    # Fill in  x and y variables
    if debug:
        print('   DEBUG main() Xsize_reg: ', Xsize_reg, " len(xarr)=", len(xarr))
    rootgrp.variables['y'][:] = yarr[:]
    rootgrp.variables['x'][:] = xarr[:]
    del yarr, xarr

    if has_subgrid:
        xarr_sub, yarr_sub = getXY(GeoTransform_subgrid, Xsize_alt, Ysize_alt, flipY=True)           # Use affine transformation to calculate cell coordinates
        if debug:
            print('   DEBUG length of yarr_sub: {l}'.format(l=len(yarr_sub)))
        rootgrp.variables['y_alt'][:] = yarr_sub[:]
        rootgrp.variables['x_alt'][:] = xarr_sub[:]
        del yarr_sub, xarr_sub
        
    # Add additional metadata to the output file
    rootgrp.GDAL_DataType = 'Generic'
    rootgrp.setncattr('WINDOW_append_vars', ','.join(append_vars))
    rootgrp.history = 'Created %s' %time.ctime()
    rootgrp.Source_Software = 'WINDOW AIR Processor %s' %PpVersion
    rootgrp.proj4 = proj4
    rootgrp.processing_notes = processing_notes_SM
    rootgrp.close()
    del rootgrp, rootgrp_in
    
    #uninstall error handler
    gdal.PopErrorHandler()
    print('Process completed in {t: 3.2f} seconds [{o}].'.format(t=(time.time()-tic), o=out_nc))