#!/usr/bin/env python3

import os, sys
import json
import glob
from optparse import OptionParser

from netCDF4 import Dataset as nc4_Dataset
from netCDF4 import Variable as nc4_Variable

from cfchecks import CFChecker, vn1_6
from cfchecks import getargs as check_getargs

WINDOW_CF_NAMES_CFG='cf_names.json'

DEBUG_P  = "   DEBUG]"
ERROR_P  = "   ERROR:"
INFO_P   = "    INFO]"
WARN_P   = "    WARN:"
ACTION_P = "  ACTION]"
FIXME_P  = "  == FIX ME =="

OPT_GIS_INPUT_Required = False

def create_parser():
    usage_str = "%prog [options] "
    parser = OptionParser(usage = usage_str)
    parser.add_option("-o", "--output-name", "--output_name", dest="out_name", default=None,
            help=" output name. must have nc extension, optional")
    parser.add_option("--out_dir", "--out-dir", dest="out_dir", default='.',
            help=" output directory - optional")
    parser.add_option("--gis_input", "--gis-input", dest="gis_input_nc", default=None,
            help=" output directory - optional")
    #Options for CFChecker
    parser.add_option('-a','--area_types', dest="areatypes", default=None)
    parser.add_option('-b','--badc', dest="badc", default=None)
    parser.add_option('-c','--coards', dest="coards", default=None)
    parser.add_option("-d", dest="debug", action="store_true", default=False,
            help=" Enable debug - optional")
    parser.add_option('-l','--uploader', dest="uploader", default=None)
    parser.add_option('-n','--noname', dest="useFileName", default=None)
    #        useFileName="no"
    parser.add_option('-s','--cf_standard_names', dest="standardname", default=None)
    parser.add_option('-v','--version', dest="cf_version", default=None)

    return parser


def correct_global_fatal(from_nc, to_nc, check_code, message):
    print('{p}   correct_global_fatal() code: {c}, {m}'.format(p=FIXME_P, c=check_code, m=message))

def correct_global_error(from_nc, to_nc, check_code, message):
    print('{p}   correct_global_error() code: {c}, {m}'.format(p=FIXME_P, c=check_code, m=message))

def correct_global_warn(from_nc, to_nc, check_code, message):
    if check_code == '2.6.1' or check_code == '(2.6.1':
        attr_key = 'Conventions'
        attr_value = vn1_6.__str__()
        to_nc.setncattr(attr_key, attr_value)
        print("{p} The global attribute {k} ({v}) is added".format(
                p=ACTION_P, k=attr_key, v=attr_value))
    else:
        print('{p}  correct_global_warn() code: {c}, message: {m}'.format(p=FIXME_P, c=check_code, m=message))

def correct_global_info(from_nc, to_nc, check_code, message):
    print('{p}  correct_global_info() code: {c}, message: {m}'.format(p=FIXME_P, c=check_code, m=message))

def correct_global_version(from_nc, to_nc, check_code, message):
    print('{p} correct_global_version() code: {c}, message: {m}'.format(p=FIXME_P, c=check_code, m=message))

def correct_variable_fatal(var_name, to_var, check_code, message, cf_name_map=None):
    print('{p}   correct_variable_fatal() {v} code: {c}, {m}'.format(
            p=FIXME_P, v=var_name, c=check_code, m=message))

def correct_variable_error(var_name, to_var, check_code, message, cf_name_map=None):
    #debug = False
    #debug = not debug
    if check_code == '3.1' and 0 == message.find('Units'):
        cf_units = cf_name_map.get_cf_units(var_name)
        if cf_units is not None:
            attr_key = 'units'
            to_var.setncattr('units', cf_units)
            print("{p} {k} attribute is changed to {v}".format(
                    p=ACTION_P, k=attr_key, v=cf_units))
    else:
        print('{p}  correct_global_error() {v} code: {c}, {m}'.format(
                p=FIXME_P, v=var_name, c=check_code, m=message))

def correct_variable_warn(var_name, to_var, check_code, message, cf_name_map):
    debug = False
    #debug = not debug
    if check_code == '3':
        nc_attrs = to_var.ncattrs()
        standard_name = cf_name_map.get_cf_standard_name(var_name)
        attr_key = 'standard_name'
        has_standard_name = attr_key in nc_attrs
        if standard_name is not None and not has_standard_name:
            to_var.setncattr(attr_key, standard_name)
            print("{p} The variable attribute {k} ({v}) is added".format(
                    p=ACTION_P, k=attr_key, v=standard_name))
            attr_key = 'units'
            if standard_name == 'fire_area' and attr_key in nc_attrs:
                nc_units = to_var.getncattr(attr_key)
                cf_units = cf_name_map.get_cf_units(var_name)
                if cf_units is not None and nc_units != cf_units:
                    to_var.setncattr(attr_key, cf_units)
                    print("{p} The variable attribute {k} is changed to {v}".format(
                            p=ACTION_P, k=attr_key, v=cf_units))
                
        long_name = cf_name_map.get_cf_long_name(var_name)
        if debug:
            print("DEBUG has_standard_name {e}, standard_name: {s}, long_name: {l}".format(
                    e=has_standard_name, s=standard_name, l=long_name))
        if long_name is None and standard_name is None:
            if debug:
                print("  DEBUG    Fix me !!!  missing  standard_name and long_name for {n}".format(n=var_name))
            long_name = var_name.lower()
        attr_key = 'long_name'
        if not has_standard_name and long_name is not None and \
                0 < len(long_name) and attr_key not in nc_attrs:
            to_var.setncattr(attr_key, long_name)
            print("{p} The variable attribute {k} ({v}) is added".format(
                    p=ACTION_P, k=attr_key, v=long_name))
    else:
        print('{p}  correct_global_warn() {v} code: {c}, {m}'.format(
                p=FIXME_P, v=var_name, c=check_code, m=message))

def correct_variable_info(var_name, to_var, check_code, message, cf_name_map=None):
    print('{p}   correct_variable_info() {v} code: {c}, {m}'.format(
            p=FIXME_P, v=var_name, c=check_code, m=message))

def correct_variable_version(var_name, to_var, check_code, message, cf_name_map=None):
    print('{p}   correct_variable_version() {v} code: {c}, {m}'.format(
            p=FIXME_P, v=var_name, c=check_code, m=message))
    

class nc_tools:

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
        
    SKIP_Attrs= [ "_FillValue"]
    @staticmethod
    def copy_variable(nc_out, from_var):
        to_var = nc_out.createVariable(from_var.name, from_var.dtype, from_var.dimensions)
        to_var[:] = from_var[:]
        for attr_name in from_var.ncattrs():
            if attr_name in nc_tools.SKIP_Attrs:
                print( " skipped attribute {v}: {a}={av}".format(
                    a=attr_name, v=from_var.name, av=from_var.getncattr(attr_name)))
                continue
            to_var.setncattr(attr_name, from_var.getncattr(attr_name))
        return to_var
        
    @staticmethod
    def create_nc_from(nc_name, from_nc):
        out_nc = nc4_Dataset(nc_name, "w")
        
        nc_tools.copy_dimensions(from_nc, out_nc)
        nc_tools.copy_global_attributes(from_nc, out_nc)
        return out_nc
        

class CF_Names:

    def __init__(self, cf_names_file=None):
        if cf_names_file is not None and not os.path.exists(cf_names_file):
            print("   WARN: The CF name configuration file () does not exist.".format(f=cf_names_file))
            cf_names_file = None
        if cf_names_file is None:
            cf_names_file = os.environ.get('CF_NAME_CONFIG', WINDOW_CF_NAMES_CFG)
        if os.path.exists(cf_names_file):
            with open(cf_names_file, encoding='utf-8') as json_file:
                self.cf_names_map = json.loads(json_file.read())
                
    def get_cf_attribute(self, key):
        attribute = None
        keys = key.split('.')
        print('key: {k}, keys: {ks}'.format(k=key, ks=keys))
        if 0 < len(keys):
            #print('   key: {k}'.format(k=keys[0]))
            #print('   key: XLAT: {v}'.format(v=self.cf_names_map.get(keys[0], None)))
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
        else:
            print("   ===== Fix me !!! WINDOW  missing {v} at json".format(v=var_name))
        
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
    
    def add_GIS_attrs(self, to_var):
        nc_dims = to_var.get_dims()
        if 1 >= len(nc_dims):
            return
        
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
        has_coordinates_attr = False
        if attr_key in nc_attr_names: 
            nc_attr = to_var.getncattr(attr_key)
            has_coordinates_attr= (nc_attr is not None)

        gis_attrs = self.gis_attrs
        if 0 < len(gis_attrs):
            for attr_key in gis_attrs.keys():
                nc_key = attr_key[len('WINDOW_'):]
                if nc_key in nc_attr_names:
                    continue
                
                if nc_key != 'esri_pe_string':
                    if has_coordinates_attr:
                        continue
                
                to_var.setncattr(nc_key, gis_attrs[attr_key])
                print("adding {k}={v}".format(k=nc_key, v=gis_attrs[attr_key]))
        
    def correct(self, from_nc_name, to_nc_name, results, categories):
        debug = False
        #debug = not debug
        global_results = results["global"]
        variables_results = results["variables"]
        if debug:
            print('===============    global_results', global_results)
            print('=============== variables_results', variables_results)
        
        from_nc = nc4_Dataset(from_nc_name, "r")
        
        to_nc   = nc_tools.create_nc_from(to_nc_name, from_nc)
        
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
        for var_name in from_nc.variables.keys():
            from_var = from_nc.variables[var_name]
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
                self.add_GIS_attrs(to_var)
        
    #def global_result_fatal(self, from_nc, results, to_nc):
    #    global_results = results["global"]
    #    variables_results = results["variables"]
    #    print("global_results", global_results)
    #    for key in variables_results.keys():
    #        variable_results = variables_results.get(key)
    #        if variable_results is None:
    #            continue
    #        print("variable_results", variable_results)

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
        offset = result.find(':')
        if 0 < offset and 10 > offset:
            check_code = result[0:offset-1].strip()
            message = result[offset+1:].strip()
            
        else:
            check_code = 'unknown'
            message = result.strip()
        
        return (check_code, message)
    def set_gis_attrs(self, gis_attrs):
        self.gis_attrs = gis_attrs

    #VARIABLE_IGNORE_CODES = {
    #      '2.1': ["Filename must have .nc suffix"],
    #    '2.6.1': ["No 'Conventions' attribute present"]

class CF_aux_tools:

    @staticmethod
    def collect_gis_attributes(gis_input_name):
        gis_attrs = {}
        gis_input_nc = nc4_Dataset(gis_input_name, "r")
        global_attrs = gis_input_nc.ncattrs()
        for global_attr in global_attrs:
            if global_attr.startswith('WINDOW_'):
                gis_attrs[global_attr] = gis_input_nc.getncattr(global_attr)
        grid_mapping_key = 'WINDOW_grid_mapping'
        if grid_mapping_key in gis_attrs.keys():
            grid_mapping_var = gis_attrs[grid_mapping_key]
            
        del gis_input_nc 
        return gis_attrs
        
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

    parser = create_parser()
    (options, args) = parser.parse_args()
    
    #Filter out arguments which are not for cfchecks
    my_args = []
    skip_next = False
    for arg in sys.argv[:]:
        if skip_next:
            skip_next = False
            continue
        if arg in ("-o", "--output-name", "--out_dir", "--out-dir", "--gis-input", "--gis-input"):
            skip_next = True
            continue
        my_args.append(arg)

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
            print('{p} The "-o" or "--output-name" options is not allowed for multiple input files.\nPlease set --out-dir option'.format(
                    p=ERROR_P))
            sys.exit(-1)
        out_files = [options.out_name]
    
    if options.gis_input_nc is None:
        print("{p} The GIS input file {g} is missing. Please set it with --gis-input".format(
                p=ERROR_P))
        if OPT_GIS_INPUT_Required:
            sys.exit(-1)
    else:
        if not os.path.exists(options.gis_input_nc):
            print("{p} The GIS input file {g} does not exist. Ignored GIS attributes".format(
                    p=WARN_P, g=options.gis_input_nc))
        else:
            gis_attrs = CF_aux_tools.collect_gis_attributes(options.gis_input_nc)
            actor.set_gis_attrs(gis_attrs)
    for in_file, out_file in zip(in_files, out_files):
        #try:
            print("in_file, out_file", in_file, out_file)
            results = checker.checker(in_file)
            actor.correct(in_file, out_file, results, checker.categories)
        #except ex:
        #    print("Processing of file %s aborted due to error" % in_file)
        #    print("ex", ex)
    
    totals = checker.get_total_counts()

    if debug:
        print("")
        print("Results dictionary:", checker.all_results)
        print("")
        print("Messages that were printed", checker.all_messages)

    errs = totals["FATAL"] + totals["ERROR"]
    if errs:
        sys.exit(errs)
    
    #warns = totals["WARN"]
    #if warns:
    #    sys.exit(-warns)

    sys.exit(0)


#--------------------------
# Main Program
#--------------------------

if __name__ == '__main__':

    main()

