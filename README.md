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
        
Command to update WRF outputs:

    python cf_convert.py -o output/wrfout_fire_d02_2018-04-17_17_00_00.nc \
    --gis-input gis_d02.nc input/wrfout_fire_d02_2018-04-17_17_00_00
    
