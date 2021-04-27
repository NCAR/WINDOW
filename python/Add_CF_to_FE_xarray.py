#!/usr/bin/env python
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# Copyright UCAR (c) 2020
# University Corporation for Atmospheric Research(UCAR)
# National Center for Atmospheric Research(NCAR)
# Research Applications Laboratory(RAL)
# P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
#
# Name:        module1
# Purpose:
# Author:      $ Kevin Sampson(ksampson)
# Created:     2021
# Licence:     <your licence>
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

'''
2021/3/19

    This script is intended to apply and alter Fast Eddy output files to contain
    roughly CF-compliant metadata. A resolution (DX,DY) must be specified, as well
    as a coordinate system origin.

Example Usage: ./Add_CF_to_FE_xarray.py ./FE_CBL.9600.nc ./FE_CBL.9600_georeferenced_wgs84_xr.nc -96.789456415999950 32.797937863000072 -96.770601578999958 32.782379042000059 86400 "seconds since 2020-10-8 15:15:42.5 -6:00"
'''

# --- Import Modules --- #
# Import Python Core Modules
import os
import sys
import time
import argparse
from math import pi

# Import Additional Modules
import numpy
import xarray as xr
from osgeo import osr
# --- End Import Modules --- #

# --- Global Variables --- #
# Dimension names for the coordinate dimensions in input file
xDim = 'xIndex'
yDim = 'yIndex'
zDim = 'zIndex'

# Coordinate variable names for the 1D x and y indices in the output file
yVar = 'y'
xVar = 'x'
zVar = 'zIndex'
tVar = 'time'

#Coordinate variable names for the input file
inXVar = "xPos"
inYVar = "yPos"

# TEMPORARY until we figure out the Fast Eddy projected coordinate system
sphere_radius = 6370000.0                                                       # Radius of sphere to use (WRF Default = 6370000.0m)
sphere_circum = 2*sphere_radius*pi
inproj = '+proj=longlat +datum=WGS84 +no_defs'                             # Proj.4 string used to define WGS84 coordinate systems. Could also use EPSG code 4326 if using ImportFromEPSG.
#inproj = '+proj=lcc +lat_1=30 +lat_2=60 +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +a={0} +b={0} +units=m +no_defs'.format(sphere_radius)

# Add a 'processing_notes' global attribute to the output file
NCnotes = ''

# DO NOT EDIT BELOW THIS LINE

# Dictionary to map old dimensions to new dimensions. This will rename the coordinate variables too
mapDims = {xDim: xVar, yDim: yVar, zDim: zVar}

# Unify all coordinate system variables to have the same name ("crs"). Ths makes it easier for WRF-Hydro output routines to identify the variable and transpose it to output files
crsVarname = True                                                               # Switch to make all coordinate system variables = "crs" instead of related to the coordinate system name
crsVar = "crs"                                                                  # Expose this as a global for other functions in other scripts to use
grid_mapping = 'crs'                                                            # Added 10/13/2017 by KMS to generalize the coordinate system variable names
CFConv = 'CF-1.5'                                                               # CF-Conventions version to place in the 'Conventions' attribute of RouteLink files
map_pro = 0                                                                     # WRF-style MAP_PROJ numbers. 0 is custom geocentric lat/lon
outNCType = 'NETCDF4'                                                           # Data model for output netCDF data ['NETCDF4', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC'] or  format=rootgrp1.data_model
# --- End Global Variables --- #

# --- Functions --- #
def add_CRS_var(sr, map_pro, CoordSysVarName, grid_mapping, PE_string, xVar, yVar, GeoTransformStr=None):
    '''
    10/13/2017 (KMS):
        This function was added to generalize the creating of a CF-compliant
        coordinate reference system variable. This was modularized in order to
        create CRS variables for both gridded and point time-series CF-netCDF
        files.
    '''
    tic1 = time.time()

    # Scalar projection variable - http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
    ##    #ds.expand_dims('')                  # Create new dimension
    ##    #ds[CoordSysVarName] = ''            # Build a new variable (scalar)
    ##    #ds[CoordSysVarName] = xr.DataArray(numpy.array(b'', dtype='|S1'))            # Build a new variable (scalar)
    ##    ds[CoordSysVarName] = xr.DataArray(numpy.array(b'', dtype='|S1'))            # Build a new variable (scalar)
    ds = xr.Dataset()
    ds[CoordSysVarName] = xr.DataArray(numpy.array(b'', dtype='|S1'))            # Build a new variable (scalar)

    # Set variable attributes
    ds[CoordSysVarName].attrs['transform_name'] = grid_mapping                                      # grid_mapping. grid_mapping_name is an alias for this
    ds[CoordSysVarName].attrs['grid_mapping_name'] = grid_mapping                                   # for CF compatibility
    ds[CoordSysVarName].attrs['esri_pe_string'] = PE_string                                         # For ArcGIS. Not required if esri_pe_string exists in the 2D variable attributes
    ds[CoordSysVarName].attrs['spatial_ref'] = PE_string                                            # For GDAl
    ds[CoordSysVarName].attrs['long_name'] = "CRS definition"                                       # Added 10/13/2017 by KMS to match GDAL format
    ds[CoordSysVarName].attrs['longitude_of_prime_meridian'] = 0.0                                  # Added 10/13/2017 by KMS to match GDAL format
    if GeoTransformStr is not None:
        ds[CoordSysVarName].attrs['GeoTransform'] = GeoTransformStr                                 # For GDAl - GeoTransform array

    # Projection specific parameters - http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
    if map_pro == 1:
        # Lambert Conformal Conic

        # Required transform variables
        ds[CoordSysVarName].attrs['_CoordinateAxes'] = '{0} {1}'.format(yVar, xVar)                 # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        ds[CoordSysVarName].attrs['_CoordinateTransformType'] = "Projection"
        ds[CoordSysVarName].attrs['standard_parallel'] = sr.GetProjParm("standard_parallel_1"), sr.GetProjParm("standard_parallel_2")     # Double
        ds[CoordSysVarName].attrs['longitude_of_central_meridian'] = sr.GetProjParm("central_meridian")     # Double. Necessary in combination with longitude_of_prime_meridian?
        ds[CoordSysVarName].attrs['latitude_of_projection_origin'] = sr.GetProjParm("latitude_of_origin")   # Double

        # Optional tansform variable attributes
        ds[CoordSysVarName].attrs['false_easting'] = sr.GetProjParm("false_easting")                # Double  Always in the units of the x and y projection coordinates
        ds[CoordSysVarName].attrs['false_northing'] = sr.GetProjParm("false_northing")              # Double  Always in the units of the x and y projection coordinates
        ds[CoordSysVarName].attrs['earth_radius'] = sphere_radius                                   # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km.
        ds[CoordSysVarName].attrs['semi_major_axis'] = sphere_radius                                # Added 10/13/2017 by KMS to match GDAL format
        ds[CoordSysVarName].attrs['inverse_flattening'] = float(0)                                  # Added 10/13/2017 by KMS to match GDAL format: Double - optional Lambert Conformal Conic parameter

    elif map_pro == 2:
        # Polar Stereographic

        # Required transform variables
        ds[CoordSysVarName].attrs['_CoordinateAxes'] = '{0} {1}'.format(yVar, xVar)                 # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        ds[CoordSysVarName].attrs['_CoordinateTransformType'] = "Projection"
        ds[CoordSysVarName].attrs['longitude_of_projection_origin'] = sr.GetProjParm("longitude_of_origin")   # Double - proj_var.straight_vertical_longitude_from_pole = ''
        ds[CoordSysVarName].attrs['latitude_of_projection_origin'] = sr.GetProjParm("latitude_of_origin")     # Double
        ds[CoordSysVarName].attrs['scale_factor_at_projection_origin'] = sr.GetProjParm("scale_factor")      # Double

        # Optional tansform variable attributes
        ds[CoordSysVarName].attrs['false_easting'] = sr.GetProjParm("false_easting")                         # Double  Always in the units of the x and y projection coordinates
        ds[CoordSysVarName].attrs['false_northing'] = sr.GetProjParm("false_northing")                       # Double  Always in the units of the x and y projection coordinates
        ds[CoordSysVarName].attrs['earth_radius'] = sphere_radius                                   # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km.
        ds[CoordSysVarName].attrs['semi_major_axis'] = sphere_radius                                # Added 10/13/2017 by KMS to match GDAL format
        ds[CoordSysVarName].attrs['inverse_flattening'] = float(0)                                  # Added 10/13/2017 by KMS to match GDAL format: Double - optional Lambert Conformal Conic parameter

    elif map_pro == 3:
        # Mercator

        # Required transform variables
        ds[CoordSysVarName].attrs['_CoordinateAxes'] = '{0} {1}'.format(yVar, xVar)                 # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        ds[CoordSysVarName].attrs['_CoordinateTransformType'] = "Projection"
        ds[CoordSysVarName].attrs['longitude_of_projection_origin'] = sr.GetProjParm("central_meridian")   # Double
        ds[CoordSysVarName].attrs['latitude_of_projection_origin'] = sr.GetProjParm("latitude_of_origin")     # Double
        ds[CoordSysVarName].attrs['standard_parallel'] = sr.GetProjParm("standard_parallel_1")                # Double
        ds[CoordSysVarName].attrs['earth_radius'] = sphere_radius                                   # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km.
        ds[CoordSysVarName].attrs['semi_major_axis'] = sphere_radius                                # Added 10/13/2017 by KMS to match GDAL format
        ds[CoordSysVarName].attrs['inverse_flattening'] = float(0)                                  # Added 10/13/2017 by KMS to match GDAL format: Double - optional Lambert Conformal Conic parameter

    elif map_pro == 6:
        # Cylindrical Equidistant or rotated pole

        #http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#appendix-grid-mappings
        # Required transform variables
        #ds[CoordSysVarName]..attrs['grid_mapping_name'] = "latitude_longitude"                      # or "rotated_latitude_longitude"

        #print('        Cylindrical Equidistant projection not supported.')
        #raise SystemExit
        pass                                                                    # No extra parameters needed for latitude_longitude

    # Added 10/13/2017 by KMS to accomodate alternate datums
    elif map_pro == 0:
        ds[CoordSysVarName].attrs['_CoordinateAxes'] = 'lat lon'
        ds[CoordSysVarName].attrs['semi_major_axis'] = sr.GetSemiMajor()
        ds[CoordSysVarName].attrs['semi_minor_axis'] =  sr.GetSemiMinor()
        ds[CoordSysVarName].attrs['inverse_flattening'] = sr.GetInvFlattening()
        pass

    # Global attributes related to CF-netCDF
    ds.attrs['Conventions'] = CFConv                                            # Maybe 1.0 is enough?
    return ds

def create_CF_NetCDF_XR(ds, srs, xVar, yVar, zVar, GTstr, ref_time, time_units_str, DX=1, DY=1, notes=''):
    """This function will create the netCDF file with CF conventions for the grid
    description. Valid output formats are 'GEOGRID', 'ROUTING_GRID', and 'POINT'.
    The output NetCDF will have the XMAP/YMAP created for the x and y variables
    and the LATITUDE and LONGITUDE variables populated from the XLAT_M and XLONG_M
    variables in the GEOGRID file or in the case of the routing grid, populated
    using the getxy function."""

    tic1 = time.time()
    print('  Creating CF-netCDF File.')

    # Build Esri WKT Projection string to store in CF netCDF file
    projEsri = srs.Clone()                                            # Copy the SRS
    projEsri.MorphToESRI()                                                      # Alter the projection to Esri's representation of a coordinate system
    PE_string = projEsri.ExportToWkt().replace("'", '"')                        # INVESTIGATE - this somehow may provide better compatability with Esri products?
    print('    Esri PE String: {0}'.format(PE_string))

    # Create coordinate variables
    var_y = ds[yVar]
    var_x = ds[xVar]
    var_z = ds[zVar]

    # Must handle difference between ProjectionCoordinateSystem and LatLonCoordinateSystem
    if srs.IsGeographic():
        if crsVarname:
            CoordSysVarName = crsVar
        else:
            CoordSysVarName = "LatLonCoordinateSystem"

        # Set variable attributes
        #var_y.attrs['standard_name'] = ''
        #var_x.attrs['standard_name'] = ''
        var_y.attrs['long_name'] = "latitude coordinate"
        var_x.attrs['long_name'] = "longitude coordinate"
        var_y.attrs['units'] = "degrees_north"
        var_x.attrs['units'] = "degrees_east"
        var_y.attrs['_CoordinateAxisType'] = "Lat"
        var_x.attrs['_CoordinateAxisType'] = "Lon"

    elif srs.IsProjected():
        if crsVarname:
            CoordSysVarName = crsVar
        else:
            CoordSysVarName = "ProjectionCoordinateSystem"
        #proj_units = sr.linearUnitName.lower()                                  # sr.projectionName wouldn't work for a GEOGCS
        proj_units = 'm'                                                        # Change made 11/3/2016 by request of NWC

        # Set variable attributes
        var_y.attrs['standard_name'] = 'projection_y_coordinate'
        var_x.attrs['standard_name'] = 'projection_x_coordinate'
        var_y.attrs['long_name'] = 'y coordinate of projection'
        var_x.attrs['long_name'] = 'x coordinate of projection'
        var_y.attrs['units'] = proj_units                                       # was 'meter', now 'm'
        var_x.attrs['units'] = proj_units                                       # was 'meter', now 'm'
        ##        var_y.attrs['_CoordinateAxisType'] = "GeoY"                             # Use GeoX and GeoY for projected coordinate systems only
        ##        var_x.attrs['._CoordinateAxisType'] = "GeoX"                            # Use GeoX and GeoY for projected coordinate systems only
        var_y.attrs['resolution'] = float(abs(DY))                              # Added 11/3/2016 by request of NWC
        var_x.attrs['resolution'] = float(DX)                                   # Added 11/3/2016 by request of NWC

        # Build coordinate reference system variable
        # (sr,GeoTransformStr) = (srs, GTstr)
        ds2 = add_CRS_var(srs, map_pro, CoordSysVarName, grid_mapping, PE_string, xVar, yVar, GeoTransformStr=GTstr)
        ds = xr.merge([ds, ds2])

    # Handle a z coordinate
    var_z.attrs['units'] = "level"
    var_z.attrs['positive'] = "down"

    # Apply metadata on a per-variable basis. This includes 'coordinates', 'grid_mapping', etc.
    for varname in list(ds.variables.keys()):
        varDims = ds[varname].dims

        # Start the list of coordinates for each variable.
        coordinates = ''

        # Allow a Z coordinate. Not used at the moment.
        ##        if zDim in varDims:
        ##            if len(coordinates) > 0:
        ##                coordinates += ' '
        ##            coordinates += '{0}'.format(zVar)

        if mapDims.get(yDim, yDim) in varDims and mapDims.get(xDim, xDim) in varDims:
            if len(coordinates) > 0:
                coordinates += ' '
            coordinates += '{0} {1}'.format(yVar, xVar)

            if srs.IsProjected():
                ds[varname].attrs['grid_mapping'] = CoordSysVarName
                #ds[varname].attrs['spatial_ref'] = PE_string                    # For GDAl
            ds[varname].attrs['coordinates'] = coordinates


    #Add CF compliant time metadata
    ds[tVar] = numpy.array([ref_time])
    var_t = ds[tVar]
    var_t.attrs['long_name'] = 'time'
    var_t.attrs['units'] = time_units_str
    var_t.attrs['calendar'] = 'gregorian' 
    # Global attributes
    ds.attrs['GDAL_DataType'] = 'Generic'
    ds.attrs['Source_Software'] = 'Custom script to test CF metadata (Kevin Sampson, NCAR)'
    ds.attrs['proj4'] = inproj
    ds.attrs['history'] = 'Created {0}'.format(time.ctime())
    ds.attrs['processing_notes'] = notes
    ds.attrs['spatial_ref'] = PE_string                                         # For GDAl
    ds.attrs['GeoTransform'] = GTstr
    print('  netCDF global attributes set after {0: 3.2f} seconds.'.format(time.time()-tic1))
    return ds

# --- End Functions --- #

# --- Main Codeblock --- #
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("in_nc",type=str,help="Input netcdf file")
    parser.add_argument("out_nc",type=str,help="Output netcdf file")
    parser.add_argument("ul_x",type=float,help="Upper left x coord")
    parser.add_argument("ul_y",type=float,help="Upper left y coord")
    parser.add_argument("lr_x",type=float,help="Lower right x coord")
    parser.add_argument("lr_y",type=float,help="Lower right y coord")
    parser.add_argument("ref_time",type=float,help="Reference time")
    parser.add_argument("time_units_str",type=str,help="UDUNITS time formatted string")
    parser.add_argument("--DX",type=float,help="x grid differential")
    parser.add_argument("--DY",type=float,help="y grid differential")
    args = parser.parse_args()

    tic = time.time()
    print('Process initiated at {0}'.format(time.ctime()))
    print('Input netCDF file: {0}'.format(args.in_nc))

    write_newfile = args.in_nc != args.out_nc

    if write_newfile:
        # Open input file on disk. Lazy loading
        ds = xr.open_dataset(args.in_nc,decode_times = False) # , decode_cf=False
    else:
        # !!! OVERWRITE exiting input file !!!
        # Load the entire dataset into memory. This allows you to overwrite the existing file
        with xr.open_dataset(args.in_nc,decode_times = False) as data:
            ds = data.load()

    DX = args.DX if args.DX else abs(args.ul_x - args.lr_x)/ds[inXVar].shape[3]
    DY = args.DY if args.DY else abs(args.ul_y - args.lr_y)/ds[inYVar].shape[2]
    GTstr = "{0} {1} 0 {2} 0 -{3}".format(args.ul_x, DX, args.ul_y, DY)
    print("GTstr : %s" % GTstr)
    # Rename the dimensions and variables using a dictionary
    ds = ds.rename(mapDims)

    # Modify the values in the renamed coordinate variables
    xMin, DX_GT, xskew, yMax, yskew, DY_GT = [float(item) for item in GTstr.split(' ')]
    for varname in  mapDims.values():
        # Transform index to coordinates using the affine geotransformation
        if varname == 'x':
            ds[varname] = xMin + 0.5*DX_GT + ds[varname][:]*DX_GT
        if varname == 'y':
            ds[varname] = yMax + 0.5*DY_GT + ds[varname][:]*DY_GT
        del varname

    # Get arrays shape
    nrows = ds[yVar].shape[0]
    ncols = ds[xVar].shape[0]

    # create the spatial reference for the input point CSV file, WGS84
    srs = osr.SpatialReference()
    srs.ImportFromProj4(inproj)

    # Run the main function to apply appropriate spatio-temporal metadata
    ds = create_CF_NetCDF_XR(ds, srs, xVar, yVar, zVar, GTstr, args.ref_time, args.time_units_str, DX=DX, DY=DY, notes=NCnotes)

    # Output file to disk
    encoding = {varname:{'_FillValue': None, 'zlib':True, 'complevel':2} for varname in list(ds.variables.keys())}
    ds.to_netcdf(args.out_nc, mode='w', format=outNCType, encoding=encoding)
    print('Output netCDF file: {0}'.format(args.out_nc))

    # Clean up
    ds.close()
    del ds, encoding, srs, nrows, ncols
    srs = None
    print('Process {0} completed in {1: 3.2f} seconds.'.format(sys.argv[0], time.time()-tic))
# --- End Main Codeblock --- #

if __name__ == '__main__':
    main()
