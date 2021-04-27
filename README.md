# WINDOW
Wide Interoperability Nexus for Data Organized Workflows

Command to generate geogrid file with GIS info

    python CF_Gridded.py geo_em.d02.nc gis_d02.subgrid.nc

        The dimentions of x and y replace the dimension west_east and south_north at WRF input.

    python CF_Gridded.py geo_em.d02.nc gis_d02.subgrid.nc subgrid_factor=4
    OR python CF_Gridded.py geo_em.d02.nc gis_d02.subgrid.nc 4

        subgrid_factor option to support west_east_subgrid and south_north_subgrid dimensions.
            note: the dimension west_east_subgrid = (west_east + 1) * subgrid_factor
                             OR west_east_subgrid = west_east * subgrid_factor
        The dimentions of x and y replace the dimension west_east and south_north at WRF input.
        The dimentions of x_alt and y_alt replace the dimension west_east_subgrid and south_north_subgrid at WRF input.

    python CF_Gridded.py geo_em.d02.nc gis_d02.subgrid.nc subgrid_factor=4 subgrid_first
    OR python CF_Gridded.py geo_em.d02.nc gis_d02.subgrid.nc 4 subgrid_first

        The subgrid_first option is ignored if subgrid_factor is not greater than 1.
        If subgrid_first is enabled, the dimentions of x and y replace the dimension west_east_subgrid
        and south_north_subgrid at WRF input. The dimentions of x_alt and y_alt replace the dimension
        west_east and south_north at WRF input.
   
Command to generate Fast Eddy files with GIS data

    usage: Add_CF_to_FE_xarray.py [-h] [--DX DX] [--DY DY] in_nc out_nc ul_x ul_y lr_x lr_y ref_time time_units_str

    positional arguments:
      in_nc           Input netcdf file
      out_nc          Output netcdf file
      ul_x            Upper left x coord
      ul_y            Upper left y coord
      lr_x            Lower right x coord
      lr_y            Lower right y coord
      ref_time        Reference time
      time_units_str  UDUNITS time formatted string

    optional arguments:
      -h, --help      show this help message and exit
      --DX DX         x grid differential
      --DY DY         y grid differential

    Example Usage:
    python Add_CF_to_FE_xarray.py ./FE_CBL.9600.nc ./FE_CBL.9600_georeferenced_wgs84_xr.nc -96.789456415999950 32.797937863000072 -96.770601578999958 32.782379042000059 86400 "seconds since 2020-10-8 15:15:42.5 -6:00"

Command to update WRF outputs:

    python cf_convert.py -o output/wrfout_fire_d02_2018-04-17_17_00_00.nc \
    --gis-input gis_d02.nc input/wrfout_fire_d02_2018-04-17_17_00_00
    

The configuration:
- The library of Unidata UDUNITS (https://www.unidata.ucar.edu/software/udunits/)
- Python packages
   - cfunits (1.9)
   - netCDF4 (1.4.2)
   - GDAL (1.11.0)
