#!/usr/bin/env python3

import os, sys
import json
import glob
from optparse import OptionParser

from cfchecks import CFChecker, vn1_6
from cfchecks import getargs as check_getargs

from window_tools import nc_tools, get_simple_logger

cf_logger = get_simple_logger()

ENV_CF_NAMES_CFG    = 'CF_NAME_CONFIG'
WINDOW_CF_NAMES_CFG = 'cf_names.json'

DEBUG_P  = "   DEBUG]"
ERROR_P  = "   ERROR:"
INFO_P   = "    INFO]"
WARN_P   = "    WARN:"
ACTION_P = "  ACTION]"
FIXME_P  = "  == FIX ME =="
CHECKME_P  = "  -- CHECK ME --"

OPT_GIS_INPUT_Required = False
OPT_add_time_dim = True

DIMENSION_MAP = {
        'x':        'west_east', 
        'y':        'south_north', 
        'x_alt':    'west_east_subgrid',
        'y_alt':    'south_north_subgrid'
    } 

ALT_DIMENSION_MAP = {
        'x_alt':    'west_east', 
        'y_alt':    'south_north',
        'x':        'west_east_subgrid',
        'y':        'south_north_subgrid'
    }
    

ARGUMENT_LIST = [ "-o", "--output-name", "--out_dir", "--out-dir",
                  "--gis-input", "--gis-input",
                  "-x","--x-dim","--x_dim", "-y","--y-dim","--y_dim"
                ]


def create_parser():
    usage_str = "%prog [options] "
    parser = OptionParser(usage = usage_str)
    parser.add_option("-o", "--output-name", "--output_name", dest="out_name", default=None,
            help=" output name. must have nc extension, optional")
    parser.add_option("--out_dir", "--out-dir", dest="out_dir", default='.',
            help=" output directory - optional")
    parser.add_option("--gis_input", "--gis-input", dest="gis_input_nc", default=None,
            help=" output directory - optional")
    parser.add_option('-x','--x_dim', '--x-dim', dest="x_dim",
            help=" dimension name for longitude" , default=None)
    parser.add_option('-y','--y_dim', '--y-dim', dest="y_dim",
            help=" dimension name for latitude" , default=None)
    #Options for CFChecker
    parser.add_option('-a','--area_types', dest="areatypes", default=None)
    parser.add_option('-b','--badc', dest="badc", default=None)
    parser.add_option('-c','--coards', dest="coards", default=None)
    parser.add_option("-d", dest="debug", action="store_true", default=False,
            help=" Enable debug - optional")
    parser.add_option("--debug-level", "--debug_level", dest="debug_level", default=0,
            help=" Set the debug level- optional")
    parser.add_option('-l','--uploader', dest="uploader", default=None)
    parser.add_option('-n','--noname', dest="useFileName", default=None)
    #        useFileName="no"
    parser.add_option('-s','--cf_standard_names', dest="standardname", default=None)
    parser.add_option('-v','--version', dest="cf_version", default=None)

    return parser


def correct_global_fatal(from_nc, to_nc, check_code, message):
    method_name = "correct_global_fatal()"
    if check_code != '2.1': # (2.1) the netcdf must has ".nc" extension.
        cf_logger.info('{p}  {n} check the global config for code {c}, {m}'.format(
                p=FIXME_P, n=method_name, c=check_code, m=message))

def correct_global_error(from_nc, to_nc, check_code, message):
    method_name = "correct_global_error()"
    cf_logger.info('{p}  {n} check the global config for code {c}\n\t{m}'.format(
            p=FIXME_P, n=method_name, c=check_code, m=message))

def correct_global_warn(from_nc, to_nc, check_code, message):
    method_name = "correct_global_warn()"
    if check_code == '2.6.1' or check_code == '(2.6.1':
        attr_key = 'Conventions'
        attr_value = vn1_6.__str__()
        to_nc.setncattr(attr_key, attr_value)
        cf_logger.info("{p} The global attribute {k} ({v}) is added".format(
                p=ACTION_P, k=attr_key, v=attr_value))
    else:
        cf_logger.info('{p}  {n} check the global config for code {c}\n\t{m}'.format(
                p=FIXME_P, n=method_name, c=check_code, m=message))

def correct_global_info(from_nc, to_nc, check_code, message):
    method_name = "correct_global_info()"
    cf_logger.info('{p}  {n} check the global config for code {c}\n\t{m}'.format(
            p=CHECKME_P, n=method_name, c=check_code, m=message))

def correct_global_version(from_nc, to_nc, check_code, message):
    method_name = "correct_global_version()"
    print('{p}  {n} check the global config for code {c}\n\t{m}'.format(
            p=FIXME_P, n=method_name, c=check_code, m=message))
    
def correct_variable_fatal(var_name, to_var, check_code, message, cf_name_map=None):
    method_name = "correct_variable_fatal()"
    cf_logger.info('{p}  {n} check the variable {v} for code: {c}\n\t{m}'.format(
            p=FIXME_P, n=method_name, v=var_name, c=check_code, m=message))

def correct_variable_error(var_name, to_var, check_code, message, cf_name_map=None):
    #debug = False
    #debug = not debug
    method_name = "correct_variable_error()"
    if check_code == '3.1' and 0 == message.find('Units'):
        cf_units = None
        if cf_name_map is not None and 0 < len(cf_name_map.cf_names_map):
            cf_units = cf_name_map.get_cf_units(var_name)
        attr_key = 'units'
        if cf_units is not None:
            to_var.setncattr('units', cf_units)
            cf_logger.log("{p} {k} attribute is changed to {v}".format(
                    p=ACTION_P, k=attr_key, v=cf_units))
        else:
            cf_logger.info("{n} {k} attribute for {v} is not available".format(
                    n=method_name, k=attr_key, v=var_name))
    else:
        print('{p}  {n} check the variable {v} for code: {c}\n\t{m}'.format(
                p=FIXME_P, n=method_name, v=var_name, c=check_code, m=message))

def correct_variable_warn(var_name, to_var, check_code, message, cf_name_map):
    debug = False
    #debug = not debug
    method_name = "correct_variable_warn()"
    if check_code == '3':
        nc_attrs = to_var.ncattrs()
        standard_name = cf_name_map.get_cf_standard_name(var_name)
        attr_key = 'standard_name'
        has_standard_name = attr_key in nc_attrs
        if standard_name is not None and not has_standard_name:
            to_var.setncattr(attr_key, standard_name)
            cf_logger.log("{p} The variable attribute {k} ({v}) is added".format(
                    p=ACTION_P, k=attr_key, v=standard_name))
            attr_key = 'units'
            if standard_name == 'fire_area' and attr_key in nc_attrs:
                #DeprecationWarning: tostring() is deprecated. Use tobytes() instead.
                nc_units = to_var.getncattr(attr_key)
                cf_units = cf_name_map.get_cf_units(var_name)
                if cf_units is not None and nc_units != cf_units:
                    to_var.setncattr(attr_key, cf_units)
                    cf_logger.log("{p} The variable attribute {k} is changed to {v}".format(
                            p=ACTION_P, k=attr_key, v=cf_units))
                
        long_name = cf_name_map.get_cf_long_name(var_name)
        cf_logger.debug(1, "has_standard_name {e}, standard_name: {s}, long_name: {l}".format(
                    e=has_standard_name, s=standard_name, l=long_name))
        if long_name is None and standard_name is None:
            cf_logger.debug(1, "Fix me !!!  missing the standard_name and long_name for {n}".format(
                    n=var_name))
            long_name = var_name.lower()
        attr_key = 'long_name'
        if not has_standard_name and long_name is not None and \
                0 < len(long_name) and attr_key not in nc_attrs:
            to_var.setncattr(attr_key, long_name)
            cf_logger.log("{p} The variable attribute {k} ({v}) is added".format(
                    p=ACTION_P, k=attr_key, v=long_name))
    else:
        cf_logger.warning('{p}  {n} check the variable {v} for code: {c}\n\t{m}'.format(
                p=FIXME_P, n=method_name, v=var_name, c=check_code, m=message))

def correct_variable_info(var_name, to_var, check_code, message, cf_name_map=None):
    method_name = "correct_variable_info()"
    if check_code == '3.1' and 0 == message.find('Units'):
        attr_key = 'units'
        is_not_corrected = True
        #if cf_name_map is not None and 0 < len(cf_name_map.cf_names_map):
        if cf_name_map is not None and 0 < len(cf_name_map.cf_names_map):
            cf_units = cf_name_map.get_cf_units(var_name)
            if cf_units is not None:
                to_var.setncattr('units', cf_units)
                is_not_corrected = False
                cf_logger.log("{p} {k} attribute is changed to {v}".format(
                        p=ACTION_P, k=attr_key, v=cf_units))
        if is_not_corrected:
            cf_logger.info("{n} {k} attribute for {v} is not available".format(
                    n=method_name, k=attr_key, v=var_name))
    else:
        cf_logger.info('{p}  {n} check the variable {v} for code: {c}\n\t{m}'.format(
                p=CHECKME_P, n=method_name, v=var_name, c=check_code, m=message))

def correct_variable_version(var_name, to_var, check_code, message, cf_name_map=None):
    method_name = "correct_variable_version()"
    cf_logger.warning('{p}  {n} check the variable {v} for code: {c}\n\t{m}'.format(
            p=FIXME_P, n=method_name, v=var_name, c=check_code, m=message))
    

class CF_Names:

    def __init__(self, cf_names_file=None):
        method_name = "CF_Names.__init__()"
        if cf_names_file is not None and not os.path.exists(cf_names_file):
            cf_logger.warning("   WARN: The CF name configuration file () does not exist.".format(
                    f=cf_names_file))
            cf_names_file = None
        if cf_names_file is None:
            cf_names_file = os.environ.get(ENV_CF_NAMES_CFG, WINDOW_CF_NAMES_CFG)
        if os.path.exists(cf_names_file):
            with open(cf_names_file, encoding='utf-8') as json_file:
                self.cf_names_map = json.loads(json_file.read())
        else:
            self.cf_names_map = {}
            cf_logger.info("{m} The configuration file '{f}' is missing which contains the standard and long names.".format(
                    m=method_name, f=WINDOW_CF_NAMES_CFG))
            cf_logger.info("\tIt should be in the current working directory or specified by using the environment variable '{e}'.".format(
                    e=ENV_CF_NAMES_CFG))
        self.missing_var_list = []
                
    def get_cf_attribute(self, key):
        attribute = None
        keys = key.split('.')
        self.logger.debug(1, 'key: {k}, keys: {ks}'.format(k=key, ks=keys))
        if 0 < len(keys):
            attribute = self.cf_names_map.get(keys[0], None)
            if attribute is not None and 1 < len(keys):
                for key in keys[1:]:
                    attribute = attribute.get(key, None)
                    if attribute is None:
                        break
        return attribute

    def get_cf_var_attribute(self, var_name, attr_key, convert_empty_to_none=True):
        attribute = self.cf_names_map.get(var_name, None)
        if attribute is not None:
            attribute = attribute.get(attr_key, None)
        elif var_name not in self.missing_var_list:
            self.missing_var_list.append(var_name)
            cf_logger.info("get_cf_var_attribute() missing {v} at json".format(v=var_name))
        
        if convert_empty_to_none and attribute is not None and 0 == len(attribute):
            attribute = None
        return attribute
      
    def get_cf_long_name(self, var_name, convert_empty_to_none=True):
        return self.get_cf_var_attribute(var_name, 'long_name')
        
    def get_cf_standard_name(self, var_name, convert_empty_to_none=True):
        return self.get_cf_var_attribute(var_name, 'standard_name')

    def get_cf_units(self, var_name, convert_empty_to_none=True):
        return self.get_cf_var_attribute(var_name, 'units')
        

class CF_Corrector:

    #VARIABLE_IGNORE_CODES = {
    #      '2.1': ["Filename must have .nc suffix"],
    #    '2.6.1': ["No 'Conventions' attribute present"]
    #}
    
    RENAME_DIMENSION = False
    
    GLOBAL_HANDLERS = {
          "FATAL": correct_global_fatal,
          "ERROR": correct_global_error,
           "WARN": correct_global_warn,
           "INFO": correct_global_info,
        "VERSION": correct_global_version
    }
    
    VARIABLE_HANDLERS = {
          "FATAL": correct_variable_fatal,
          "ERROR": correct_variable_error,
           "WARN": correct_variable_warn,
           "INFO": correct_variable_info,
        "VERSION": correct_variable_version
    }
    
    def __init__(self):
        self.gis_attrs = {}
        self.cf_name_map = CF_Names()
        self.gis_dims = []
        self.gis_vars  = []
        self.subgrid_first = False
        self.new_dims = DIMENSION_MAP
        self.logger = get_simple_logger()
    
    def add_gis_attrs(self, to_var):
        method_name = 'CF_Corrector.add_gis_attrs()'
        debug_level = 3
        #nc_dims = to_var.get_dims()
        nc_dims = to_var.dimensions
        if 1 >= len(nc_dims):
            return
        
        has_sub_dim = False
        has_time_dim = False
        
        alt_dim_names = []
        new_dim_name = self.new_dims.get('x_alt', None)
        if new_dim_name is not None:
            alt_dim_names.append(new_dim_name)
        new_dim_name = self.new_dims.get('y_alt', None)
        if new_dim_name is not None:
            alt_dim_names.append(new_dim_name)
        
        #for dim in nc_dims:
        #    dim_name = dim.name
        for dim_name in nc_dims: 
            if 'time' == dim_name.lower() or 'xtime' == dim_name.lower():
                has_time_dim = True
            elif dim_name.endswith('_alt') or dim_name in alt_dim_names:
                has_sub_dim = True
        
        var_name = to_var.name
        nc_attr_names = to_var.ncattrs()
        
        attr_key = 'long_name'
        if attr_key in nc_attr_names: 
            nc_attr = to_var.getncattr(attr_key)
            if nc_attr is not None:
                if nc_attr == 'times':
                    return
        
        attr_key = 'standard_name'
        if attr_key in nc_attr_names: 
            nc_attr = to_var.getncattr(attr_key)
            if nc_attr is not None:
                if nc_attr == 'time':
                    return
        
        attr_key = 'coordinates'
        #has_coordinates_attr = False
        if attr_key in nc_attr_names: 
            nc_attr = to_var.getncattr(attr_key)
            has_coordinates_attr = (nc_attr is not None)
            if has_coordinates_attr:
                self.logger.debug(debug_level, "{m} {k}={v} for {vn}".format(
                        m=method_name, k=attr_key, v=nc_attr,vn=var_name))
                new_nc_attr = nc_attr
                #if "XLONG XLAT XTIME" == nc_attr and 'time' == nc_dims[0].name.lower():
                if "XLONG XLAT XTIME" == nc_attr:
                    #nc_attr = "XTIME XLAT XLONG"
                    new_nc_attr = "XTIME XLAT XLONG"
                elif "XLONG XLAT" == nc_attr:
                    #nc_attr = "XTIME XLAT XLONG"
                    new_nc_attr = "XLAT XLONG"
                if has_sub_dim:
                    new_nc_attr = new_nc_attr.replace('XLAT', 'y_alt').replace('XLONG', 'x_alt')
                else:
                    new_nc_attr = new_nc_attr.replace('XLAT', 'y').replace('XLONG', 'x')
                
                if new_nc_attr != nc_attr:
                    to_var.setncattr(attr_key, new_nc_attr)
                    self.logger.debug(debug_level, "{m} Updated {k}={v} for {vn}".format(
                            m=method_name, k=attr_key, v=new_nc_attr,vn=var_name))

        gis_attrs = self.gis_attrs
        if 0 < len(gis_attrs):
            for attr_key in gis_attrs.keys():
                if 'WINDOW_append_vars' == attr_key:
                    continue
                
                attr_value = gis_attrs[attr_key]
                nc_raw_key = attr_key[len('WINDOW_'):]
                is_alt_key = False
                nc_key = nc_raw_key
                if 'grid_mapping_alt' == nc_raw_key:
                    is_alt_key = True
                    nc_key = 'grid_mapping'
                elif 'CoordinateSystems_alt' == nc_raw_key:
                    is_alt_key = True
                    nc_key = 'CoordinateSystems'
                
                if nc_key in nc_attr_names:
                    self.logger.debug(debug_level, '{m} Exist attribute {k} for {vn} already'.format(
                            m=method_name, k=nc_key, vn=var_name))
                    continue
                
                self.logger.debug(debug_level, '{m} gis_attr: {ko} => {k} = {v} for {vn}'.format(
                        m=method_name, ko=nc_raw_key, k=nc_key, v=attr_value, vn=var_name))
                if ('grid_mapping' == nc_key or 'CoordinateSystems' == nc_key):
                    if is_alt_key != has_sub_dim:
                        self.logger.debug(debug_level,
                                '{m} Ignored attribute {k} because is_alt_key={a} vs. has_sub_dim={s} for {vn}'.format(
                                m=method_name, k=nc_raw_key, a=is_alt_key, s=has_sub_dim, vn=var_name))
                    
                    if has_time_dim and OPT_add_time_dim:
                        attr_value = attr_value + '_t'
                #elif nc_key != 'esri_pe_string':
                #    if has_coordinates_attr:
                #        continue
                
                to_var.setncattr(nc_key, attr_value)
                self.logger.debug(debug_level,
                        "{m} adding {k}={v} for {vn}".format(
                        m=method_name, k=nc_key, v=attr_value,vn=var_name))
        
#     def add_gis_dimensions(self, to_nc):
#         method_name = "add_gis_dimensions()"
#         debug = False
#         debug = not debug
# 
#         global_dims = to_nc.dimensions
#         if self.gis_dim_x.name not in global_dims.keys():
#             to_nc.createDimension(self.gis_dim_x.name, self.gis_dim_x.size)
#         if self.gis_dim_y.name not in global_dims.keys():
#             to_nc.createDimension(self.gis_dim_x.name, self.gis_dim_x.size)
        
    def add_gis_vars(self, to_nc, new_dims={}):
        method_name = "add_gis_vars()"
        debug_level = 3
        if self.logger.is_log_enabled(debug_level):
            self.logger.debug(3, "{m} is called, global_dims:".format(m=method_name))
            for a_dim in to_nc.dimensions:
                self.logger.debug(3, "\t{d}".format(d=a_dim))
        #if self.gis_dim_x is not None and self.gis_dim_y is not None and 0 < len(self.gis_vars):
        if 0 < len(self.gis_dims) and 0 < len(self.gis_vars):
            for gis_var in self.gis_vars:
                nc_tools.copy_variable(to_nc, gis_var, new_dims)
        
    def correct(self, from_nc_name, to_nc_name, results, categories):
        method_name = "CF_Corrector.correct()"
        global_results = results["global"]
        variables_results = results["variables"]
        self.logger.debug(1, '===============    global_results {v}'.format(v=global_results))
        self.logger.debug(1, '=============== variables_results {v}'.format(v=variables_results))
        
        from_nc = nc_tools.open_nc_Dataset(from_nc_name)
        to_nc   = nc_tools.create_nc_from(to_nc_name, from_nc, self.gis_dims)
        self.logger.debug(1,
                ' correct()  new global dimensions: {v}'.format(v=to_nc.dimensions))
        
        for category in categories:
            if category == 'VERSION':
                continue
            category_handler = CF_Corrector.GLOBAL_HANDLERS[category]
            category_results = global_results[category]
            if 0 < len(category_results):
                self.correct_global(category_handler, from_nc, to_nc,
                        category_results)
        
        gis_attrs = self.gis_attrs
        results_var_names = variables_results.keys()
        self.logger.debug(1, 'correct() dimensions for variable: {v}'.format(v=self.new_dims))
            
        if 0 < len(self.gis_dims):
            found_error = False
            for gis_dim in self.gis_dims:
                data_dim_name = self.new_dims.get(gis_dim.name)
                if data_dim_name is not None and data_dim_name in to_nc.dimensions:
                    gis_dim_size = gis_dim.size
                    data_dim_size = to_nc.dimensions[data_dim_name].size
                    if gis_dim_size != data_dim_size:
                        found_error = True
                        cf_logger.error('{m} The dimension mismatch: {d1}: {s1} != {d2}: {s2}'.format(
                                m=method_name, d1=data_dim_name, s1=data_dim_size,
                                d2=gis_dim.name, s2=gis_dim_size))
            if found_error:
                cf_logger.error('{m} Did not convert because of different dimensions (could be a different domain)'.format(
                        m=method_name))
                return
            
        for var_name in from_nc.variables.keys():
            from_var = from_nc.variables[var_name]
            if CF_Corrector.RENAME_DIMENSION and 0 < len(self.gis_dims):
                to_var = nc_tools.copy_variable(to_nc, from_var, self.new_dims)
            else:
                to_var = nc_tools.copy_variable(to_nc, from_var)
                
            if var_name in results_var_names:
                cf_results = variables_results.get(var_name)
                if cf_results is not None:
                    for category in categories:
                        if category == 'VERSION':
                            continue
                        category_results = cf_results.get(category)
                        if category_results is not None and 0 < len(category_results):
                            category_handler = CF_Corrector.VARIABLE_HANDLERS[category]
                            self.correct_variable(category_handler, var_name,
                                    to_var, category_results, self.cf_name_map)
            
            if 0 < len(gis_attrs):
                self.add_gis_attrs(to_var)
        
        self.add_gis_vars(to_nc, self.new_dims)
        
    #def global_result_fatal(self, from_nc, results, to_nc):
    #    global_results = results["global"]
    #    variables_results = results["variables"]
    #    cf_logger.log("global_results", global_results)
    #    for key in variables_results.keys():
    #        variable_results = variables_results.get(key)
    #        if variable_results is None:
    #            continue
    #        cf_logger.log("variable_results", variable_results)

    def correct_global(self, category_handler, from_nc, to_nc, results):
        for result in results:
            (check_code, message) = self.separate_error_code(result)
            category_handler(from_nc, to_nc, check_code, message)
        
    def correct_variable(self, category_handler, var_name,
                         to_var, results, cf_name_map):
        for result in results:
            (check_code, message) = self.separate_error_code(result)
            category_handler(var_name, to_var, check_code, message, cf_name_map)

    def separate_error_code(self, result):
        method_name = "separate_error_code()"
        offset = result.find(':')
        if 0 < offset and 10 > offset:
            start_offset = 1 if '(' == result[0] else 0
            end_offset = (offset-2) if ')' == result[offset] else (offset-1)
            check_code = result[start_offset:end_offset].strip()
            message = result[offset+1:].strip()
            
        else:
            check_code = 'unknown'
            message = result.strip()
        
        self.logger.debug(1, "{n} check_code: '{c}' message: [{m}]".format(
                n=method_name, c=check_code, m=message))
        return (check_code, message)
    
    def set_gis_data(self, gis_data):
        self.gis_attrs = gis_data[0]
        self.gis_dims  = gis_data[1]
        self.gis_vars  = gis_data[2]
        
        dim_x = None
        dim_x_alt = None
        for dim in self.gis_dims:
            dim_name = dim.name 
            if dim_name == 'x':
                dim_x = dim
            elif dim_name == 'x_alt':
                dim_x_alt = dim

        if dim_x is not None and dim_x_alt is not None:
            self.subgrid_first = dim_x.size > dim_x_alt.size
        
        self.new_dims = ALT_DIMENSION_MAP if self.subgrid_first else DIMENSION_MAP

    #VARIABLE_IGNORE_CODES = {
    #      '2.1': ["Filename must have .nc suffix"],
    #    '2.6.1': ["No 'Conventions' attribute present"]

class CF_aux_tools:

    @staticmethod
    def collect_gis_data(gis_input_name):
        method_name = 'collect_gis_data()'
        debug = False
        #debug = not debug

        #self.logger.debug("{m} is called".format(m=method_name))
        gis_attrs = {}
        var_names = []
        gis_input_nc = nc_tools.open_nc_Dataset(gis_input_name)
        global_attrs = gis_input_nc.ncattrs()
        for global_attr in global_attrs:
            if global_attr.startswith('WINDOW_'):
                gis_attrs[global_attr] = gis_input_nc.getncattr(global_attr)
                cf_logger.debug(1, "{m} attr: {k} = {v}".format(
                        m=method_name, k=global_attr, v=gis_attrs[global_attr]))
                if global_attr == 'WINDOW_append_vars':
                    var_names = gis_attrs[global_attr].split(',')
        #grid_mapping_key = 'WINDOW_grid_mapping'
        #if grid_mapping_key in gis_attrs.keys():
        #    grid_mapping_var = gis_attrs[grid_mapping_key]
        
        if 0 == len(var_names):
            var_names = ['x', 'y', 'crs', 'crs_t', 'LATITUDE', 'LONGITUDE' ]

        dims = []
        dim_keys = []
        variables = []
        for var_name in var_names:
            var = gis_input_nc.variables.get(var_name, None)
            if var is not None:
                variables.append(var)
                #if self.logger.debug(1, "{m} dimensions: {d}".format(m=method_name, d=var.dimensions))
                for dim_name in var.dimensions:
                    if not dim_name in dim_keys:
                        dim_keys.append(dim_name)
                        dims.append(gis_input_nc.dimensions[dim_name])
            
        del gis_input_nc
        return gis_attrs, dims, variables
        
    @staticmethod
    def make_output_names(in_files, output_dir):
        out_files = []
        for in_file in in_files:
            output_name = os.path.basename(in_file)
            if not output_name.endswith(".nc"):
                output_name = output_name + ".nc"
            if output_dir is not None:
                output_name = os.path.join(output_dir, output_name)
            out_files.append(output_name)
        return out_files
        
    @staticmethod
    def make_output_names_from_input_dir(in_dir, output_dir):
        in_files = glob.glob(os.path.join(in_dir, '*'))
        out_files = CF_aux_tools.make_output_names(in_files, output_dir)
        return out_files
    
def main():

    method_name = "main()"
    parser = create_parser()
    options_args = parser.parse_args()
    options = options_args[0]
    
    #args = options_args[1]
    cf_logger.set_log_level(options.debug_level)

    
    #Filter out arguments which are not for cfchecks
    my_args = []
    skip_next = False
    for arg in sys.argv[:]:
        if skip_next:
            skip_next = False
            continue
        if arg not in ("-d"):
            if arg in ARGUMENT_LIST or arg.split('=')[0] in ARGUMENT_LIST:
                skip_next = True
                continue
            my_args.append(arg)
    
    if options.x_dim is not None:
        DIMENSION_MAP['x'] = options.x_dim
    if options.y_dim is not None:
        DIMENSION_MAP['y'] = options.y_dim
    
    (badc,coards,uploader,useFileName,standardName,areaTypes,
            version,in_files,debug)=check_getargs(my_args)

    actor = CF_Corrector()
    checker = CFChecker(uploader=uploader, useFileName=useFileName,
                        badc=badc, coards=coards, cfStandardNamesXML=standardName,
                        cfAreaTypesXML=areaTypes, version=version, debug=debug)
    if options.out_name is None:
        out_files = CF_aux_tools.make_output_names(in_files, options.out_dir)
    else:
        if 1 < len(in_files):
            cf_logger.error('The "-o" or "--output-name" options is not allowed for multiple input files.\nPLease set --out-dir option.')
            print('in_files', in_files)
            sys.exit(-1)
        out_nc_name = options.out_name
        if not out_nc_name.endswith(".nc"):
            out_nc_name += ".nc"
        out_files = [out_nc_name]
    
    if options.gis_input_nc is None:
        cf_logger.warning("The GIS input file is missing. Please set it with --gis-input.")
        if OPT_GIS_INPUT_Required:
            sys.exit(-1)
    else:
        if not os.path.exists(options.gis_input_nc):
            cf_logger.warning("The GIS input file {g} does not exist. Ignored GIS attributes.".format(
                    g=options.gis_input_nc))
        else:
            gis_data = CF_aux_tools.collect_gis_data(options.gis_input_nc)
            actor.set_gis_data(gis_data)
    
    if options.debug:
        nc_tools.set_debug_option(options.debug)
        
    for in_file, out_file in zip(in_files, out_files):
        #try:
            cf_logger.info("{m} in_file: {i}, out_file: {o}\n".format(
                    m=method_name, i=in_file, o=out_file))
            results = checker.checker(in_file)
            cf_logger.log("\n==================")
            cf_logger.log("Start correcting {m} ...".format(m=method_name))
            cf_logger.log("==================")
            actor.correct(in_file, out_file, results, checker.categories)
        #except ex:
        #    print("Processing of file %s aborted due to error" % in_file)
        #    print("ex", ex)
    
    totals = checker.get_total_counts()

    debug_level = 5
    if cf_logger.is_log_enabled(debug_level):
        cf_logger.debug(debug_level, "")
        cf_logger.debug(debug_level, "Results dictionary: {v}".format(v=checker.all_results))
        cf_logger.debug(debug_level, "")
        cf_logger.debug(debug_level, "Messages that were printed {v}".format(v=checker.all_messages))
        cf_logger.log("")

    err_count = totals["FATAL"] + totals["ERROR"]
    if err_count:
        cf_logger.log("Done with detecting {e} critical errors".format(e=err_count))
        sys.exit(err_count)
    
    #warns = totals["WARN"]
    #if warns:
    #    sys.exit(-warns)

    cf_logger.log("Done")
    sys.exit(0)


#--------------------------
# Main Program
#--------------------------

if __name__ == '__main__':

    main()

