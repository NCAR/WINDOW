'''
Created on May 24, 2019

@author: hsoh
'''

from netCDF4 import Dataset as nc4_Dataset
#from netCDF4 import Variable as nc4_Variable

class nc_tools:
    debug = False
    
    @staticmethod
    def copy_dimensions(from_nc, to_nc):
        to_nc_dims = {}
        for dim_name in from_nc.dimensions:
            nc_dim = from_nc.dimensions[dim_name]
            #print('dim', dim_name, from_nc.dimensions[dim_name])
            to_nc_dims[dim_name] = to_nc.createDimension(dim_name, len(nc_dim))
            
        #print('to_nc_dims', to_nc_dims)
        return to_nc_dims
        
    @staticmethod
    def copy_global_attributes(from_nc, to_nc):
        for attr_name in from_nc.ncattrs():
            #print('\t%s:' % attr_name, repr(from_nc.getncattr(attr_name)))
            to_nc.setncattr(attr_name, from_nc.getncattr(attr_name))
        
    SKIP_Attrs = [ '_FillValue']
    @staticmethod
    def copy_variable(nc_out, from_var, new_dims={}, var_name=None):
        cf_logger = get_simple_logger()
        if var_name is None:
            var_name = from_var.name
        #dims = [ dim if new_dims.get(dim,None) is None else new_dims[dim] for dim in from_var.dimensions]
        from_dims = [ dim for dim in from_var.dimensions ]
        if 0 < len(new_dims.keys()):
            #dims = [ dim if ( new_dims.get(dim,None) is None
            #                or nc_out.dimensions.get(new_dims.get(dim, "not_named"), None) is None )
            #             else new_dims[dim] for dim in from_var.dimensions ]
            dims = []
            for dim in from_dims:
                to_dim = dim
                for dim_key, dim_name in new_dims.items():
                    if dim_name == dim:
                        to_dim = dim_key
                        break
                dims.append(to_dim)
            
        else:
            dims = from_dims;
        cf_logger.debug(3,' copy_variable() var: {v} dims: {d} from {o}'.format(
                v=from_var.name, d=dims, o=from_var.dimensions))
            #print('   DEBUG copy_variable() new_dims: {nd} from_dims: {fd}'.format(nd=new_dims, fd=from_dims))
        to_var = nc_out.createVariable(var_name, from_var.dtype, dims)
        if 0 < len(from_dims):
            to_var[:] = from_var[:]
        for attr_name in from_var.ncattrs():
            if attr_name in nc_tools.SKIP_Attrs:
                cf_logger.log( " skipped attribute {v}: {a}={av}".format(
                    a=attr_name, v=from_var.name, av=from_var.getncattr(attr_name)))
                continue
            to_var.setncattr(attr_name, from_var.getncattr(attr_name))
        if 'XTIME' == var_name or 'XTIME' == from_var.name:
            attr_name = '_CoordinateAxisType'
            if not attr_name in from_var.ncattrs():
                to_var.setncattr(attr_name, 'Time')
                
        return to_var
        
    @staticmethod
    def create_nc_Dataset(nc_name, nc_format='NETCDF4'):
        return nc4_Dataset(nc_name, "w", format=nc_format)
        
    @staticmethod
    def create_nc_from(nc_name, from_nc, extra_dims=None):
        out_nc = nc4_Dataset(nc_name, "w")
        
        nc_tools.copy_dimensions(from_nc, out_nc)
        if extra_dims is not None:
            global_dims = out_nc.dimensions.keys()
            for dim in extra_dims:
                if dim.name not in global_dims:
                    out_nc.createDimension(dim.name, dim.size)
        nc_tools.copy_global_attributes(from_nc, out_nc)
        return out_nc
        
    @staticmethod
    def find_lat_long_variables(nc):
        method_name = "find_lat_long_variables(nc)"
        cf_logger = get_simple_logger()
        variables = []
        attr_name = 'units'
        for var_name in nc.variables.keys():
            var = nc.variables[var_name]
            if attr_name in var.ncattrs():
                units_str = var.getncattr(attr_name)
                if units_str.lower().startswith('degree'):
                    cf_logger.debug(1, "{m} {v}.units: {u}".format(
                            m=method_name, v=var_name, u=units_str))
                    variables.append(var)
        return variables
        
    @staticmethod
    def open_nc_Dataset(nc_name):
        return nc4_Dataset(nc_name, "r")

    @staticmethod
    def set_debug_option(debug):
        nc_tools.debug = debug


class SimpleLogger(object):
    
    singleton = None
    
    def __init__(self, log_level=1):
        self.set_log_level(log_level)
        
    def is_log_enabled(self, log_level):
        return self.log_level >= log_level

    def set_log_level(self, log_level):
        self.log_level = log_level

    def debug(self, log_level, message):
        if self.is_log_enabled(log_level):
            self.log(" Debug({l}) {m}".format(l=log_level, m=message))
    
    def log(self, message):
        print(message)
    
    def error(self, message):
        self.log(" == ERROR ==")
        self.log(" == ERROR == {m}".format(m=message))
        self.log(" == ERROR ==")
    
    def info(self, message):
        self.log("    INFO] {m}".format(m=message))

    def warning(self, message):
        self.log("  -- WARN -- {m}".format(m=message))

def get_simple_logger():
    if SimpleLogger.singleton is None:
        SimpleLogger.singleton = SimpleLogger()
    return SimpleLogger.singleton