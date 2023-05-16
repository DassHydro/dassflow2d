#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
    Tools to generate python, yacs scripts from data

    Author(s): Fabrice Zaoui, Yoann Audouin, Cedric Goeury

    Copyright EDF 2016
"""
from os import linesep, path, getcwd, sep, remove, rename
import glob
from utils.polygon import import_poly_from_file
from utils.exceptions import TelemacException

# Template for the Python script of a study
STUDY_TEMPLATE =\
    """
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    # Class {module} import
    from telapy.api.t2d import {module}
    from mpi4py import MPI

    TMP_DIR_NAME = 'tmp_'

    def create_tmp_folder(tmp_path,data_path):
        "Creating temporary folder and copying input files folder"
        import string, random
        from os import sep, mkdir, walk
        from shutil import copy, copytree
        #
        liste_char = string.ascii_letters+string.digits
        compteur = 7
        random_string = ""
        for i in range(compteur):
            random_string += liste_char[random.randint(0,len(liste_char)-1)]
        tmp_dir = tmp_path + sep + TMP_DIR_NAME + str(random_string)
        mkdir(tmp_dir)

        # Copying all that is in the case dir into temporary folder
        _, subdirs, files = next(walk(data_path))
        for fle in files:
            copy(data_path+sep+fle, tmp_dir)
        for subdir in subdirs:
            copytree(data_path+sep+subdir, tmp_dir+sep+subdir)

        return tmp_dir, random_string

    def delete_tmp_folder(tmp_path, res_path, random_string, res_file):
        "Removing temporary folder and copying back result file"
        from os import path
        from shutil import copy, rmtree
        #
        copy(path.join(tmp_path, 'input.txt'),
             path.join(res_path, 'input_'+random_string+'.txt'))
        copy(path.join(tmp_path, 'output.txt'),
             path.join(res_path, 'output_'+random_string+'.txt'))
        root, ext = path.splitext(res_file)
        copy(path.join(tmp_path, res_file),
             path.join(res_path, root+'_'+random_string+ext))

        rmtree(tmp_path)

    def dump_input_data({inputs}):
        "Dumping input data into a input.txt file"
        input_file = 'input.txt'
        with open(input_file,'w') as f:
            f.write('#{input_names}\\n'.replace(',', ' '))
            for iput in [{input_names}]:
                f.write(str(iput)+' ')

    def dump_output_data({outputs}):
        "Dumping output data into a input.txt file"
        "Creating temporary folder and copying input files folder"
        output_file = 'output.txt'
        with open(output_file,'w') as f:
            f.write('#{outputs}\\n'.replace(',', ' '))
            for output in [{outputs}]:
                f.write(str(output)+' ')

    def study_function({inputs}):
        from os import path, sep, chdir
        comm = MPI.COMM_WORLD
        steering_file = path.basename('{steering_file}')
        case_dir = path.dirname('{steering_file}')
        if '{user_fortran}' == 'None':
            user_fortran = None
        else:
            user_fortran = '{user_fortran}'.replace(case_dir+sep, '')
        # Creating a temporary folder to run the study
        tmp_dir, random_string = create_tmp_folder('{working_dir}', case_dir)
        # Moving to the study folder
        chdir(tmp_dir)
        dump_input_data({input_names})
        # Creation of the instance {module}
        my_case = {module}(steering_file, user_fortran=user_fortran, comm=comm)
        {run_config_getset}
        # Reading the steering file informations
        my_case.set_case()
        {run_read_case_getset}
        # Initalization
        my_case.init_state_default()
        {run_init_getset}
        # Run all time steps
        my_case.run_all_time_steps()
        {run_timestep_getset}
        # Ending the run
        my_case.finalize()
        {run_finalize_getset}
        # Instance delete
        del(my_case)
        dump_output_data({outputs})
        # Creating a temporary folder to run the study
        delete_tmp_folder(tmp_dir, '{res_path}', random_string, res_file)

        return {outputs}

    {outputs} = study_function({inputs})
    """
NO_POSITION = -1
RUN_SET_CONFIG_POS = 0
RUN_READ_CASE_POS = 1
RUN_ALLOCATION_POS = 2
RUN_INIT_POS = 3
RUN_TIMESTEP_POS = 4
RUN_FINALIZE_POS = 5

VARINFO = {}


def jdc_to_dict(jdc, command_list):
    """
    This tricky function transforms a JdC with a single command into a
    dictionary that can be used more easily from a Python context (thanks to
    M. Courtois and G. Boulant).
    """
    context = {}
    for command in command_list:
        context[command] = _args_to_dict
    exec("parameters = " + jdc.strip() in context)
    return context['parameters']


def _args_to_dict(**kwargs):
    """
    Unreferrencing
    """
    return kwargs


def get_jdc_dict_var_as_tuple(jdc_dict, varname):
    """
    Convert jdc to tuple
    """
    if varname not in jdc_dict:
        return tuple()
    if not isinstance(jdc_dict[varname], tuple):
        return (jdc_dict[varname],)

    return jdc_dict[varname]


def generate_set(input_variable):
    """
    Will generate a string containing a call to the set of the api

    @param input_variable Eficas input_variable rubrique

    @returns a string
    """

    string = ""

    # Getting name of python variable
    name = input_variable['NAME']
    # Getting name of the variable in the API
    var_name = input_variable['VAR_INFO']['VAR_NAME']

    # Checking what type of zone definition we have
    if 'INDEX' in input_variable['VAR_INFO']['ZONE_DEF']:
        i, j, k = input_variable['VAR_INFO']['ZONE_DEF']['INDEX']
        string = "my_case.set('%s', %s, i=%i, j=%i, k=%i)" % (var_name,
                                                              name, i,
                                                              j, k)
    elif 'RANGE' in input_variable['VAR_INFO']['ZONE_DEF']:
        irange = input_variable['VAR_INFO']['ZONE_DEF']['RANGE']
        string = "my_case.set_on_range('%s', %s, '%s')" % \
                 (var_name, name, irange)
    elif 'POLYGON_FILE' in input_variable['VAR_INFO']['ZONE_DEF']:
        poly_file = \
            input_variable['VAR_INFO']['ZONE_DEF']['POLYGON_FILE']['FILE_NAME']
        poly_sep = \
            input_variable['VAR_INFO']['ZONE_DEF']['POLYGON_FILE']['SEPARATOR']
        poly = import_poly_from_file(poly_file, poly_sep)
        string = "my_case.set_on_polygon('%s', %s, %s)" % (var_name,
                                                           name, str(poly))
    elif 'POLYGON' in input_variable['VAR_INFO']['ZONE_DEF']:
        poly = input_variable['VAR_INFO']['ZONE_DEF']['POLYGON']
        string = "my_case.set_on_polygon('%s', %s, %s)" % (var_name,
                                                           name, str(poly))
    else:
        raise TelemacException("Missing Zone definition")

    return string


def generate_get(output_variable):
    """
    Will generate a string containing a call to the get of the api

    @param output_variable Eficas output_variable rubrique

    @returns a string
    """

    string = ""

    # Getting name of python variable
    name = output_variable['NAME']
    # Getting name of the variable in the API
    var_name = output_variable['VAR_INFO']['VAR_NAME']

    # Checking what type of zone definition we have
    if 'INDEX' in output_variable['VAR_INFO']['ZONE_DEF']:
        i, j, k = output_variable['VAR_INFO']['ZONE_DEF']['INDEX']
        string = "%s = my_case.get('%s', i=%i, j=%i, k=%i)" % (name,
                                                               var_name,
                                                               i, j, k)
    else:
        raise TelemacException("Missin Zone definition")

    return string


def handle_variables(jdc):
    """
    Extract all the informations for the input and ouput variables

    @param jdc (dict) eficas dictionary

    @returns list of input names, list if input types, list of input values
             list of output names, list of output types, getset array
    """

    input_names = []
    input_types = []
    input_val = []
    output_names = []
    output_types = []
    getset = ["", "", "", "", "", "", ""]

    input_vars = get_jdc_dict_var_as_tuple(jdc, "INPUT_VARIABLE")
    output_vars = get_jdc_dict_var_as_tuple(jdc, "OUTPUT_VARIABLE")

    input_names = [var["NAME"].strip() for var in input_vars]
    output_names = [var["NAME"].strip() for var in output_vars]

    for var in input_vars:
        varname = var['VAR_INFO']['VAR_NAME'].rstrip()
        pos = VARINFO[varname]['set_pos']
        vartype = VARINFO[varname]['type']
        if pos == NO_POSITION:
            raise TelemacException(["Unknown position for "+var])
        getset[pos] += generate_set(var) + linesep + ' '*4
        input_val.append(var['VAR_INFO']['DEFAULT_VALUE'])
        input_types.append(vartype)
    for var in output_vars:
        varname = var['VAR_INFO']['VAR_NAME'].rstrip()
        pos = VARINFO[varname]['get_pos']
        vartype = VARINFO[varname]['type']
        if pos == NO_POSITION:
            raise TelemacException(["Unknown position for "+var])
        output_types.append(vartype)
        getset[pos] += generate_get(var) + linesep + ' '*4

    return input_names, input_types, input_val,\
        output_names, output_types, getset


def build_var_info(jdc):
    """
    Build the VARINFO dictionary

    @param jdc The eficas dictionary
    """
    from telapy.api.t2d import Telemac2d

    with open('dummy.cas', 'w') as f:
        f.write('/Dummy steering file')

    t2d = Telemac2d('dummy.cas')

    VARINFO.update(t2d.generate_var_info())

    del t2d


def get_result_keyword(module):
    """
    Get the keyword for the result file

    @param module Module of TELEMAC-MACARET
    """
    if module == 'Telemac2d':
        key_name = 'MODEL.RESULTFILE'
    else:
        raise TelemacException(
            "module "+module+" is not handled "
            "by generate_study yet!")

    return key_name


def get_result_file(module):
    """
    Generate the Api command to set the result file name

    @param module Module of TELEMAC-MACARET
    """

    key_name = get_result_keyword(module)

    var = {'NAME': 'res_file',
           'VAR_INFO': {'VAR_NAME': key_name,
                        'ZONE_DEF': {'INDEX': (0, 0, 0)}
                        }
           }
    return generate_get(var)


def set_result_file(module, res_file):
    """
    Generate the Api command to set the result file name

    @param module Module of TELEMAC-MACARET
    @param res_file Name of the new result file
    """

    key_name = get_result_keyword(module)

    var = {'NAME': 'tmpname',
           'VAR_INFO': {'VAR_NAME': key_name,
                        'ZONE_DEF': {'INDEX': (0, 0, 0)}
                        }
           }
    return generate_set(var).replace('tmpname', "'{}'".format(res_file))


def generate_study_script(jdc):
    """
    Builds a study function from an Eficas dictionary

    @param jdc The eficas dictionary

    @returns A string containing the python script of the study
    """

    build_var_info(jdc)

    steering_file = jdc['STEERING_FILE']

    user_fortran = jdc.get('USER_FORTRAN', 'None')
    module = 'Telemac2d'

    input_names, _, input_val, output_names, _, getset = handle_variables(jdc)

    inputs = ["%s=%s" %
              input_data for input_data in zip(input_names, input_val)]

    # Adding result file if in the catalog
    res_folder = jdc.get('RESULT_DIRECTORY', path.dirname(steering_file))

    if 'RESULTS_FILE_NAME' in jdc:
        getset[RUN_READ_CASE_POS] += set_result_file(
            module,
            jdc['RESULTS_FILE_NAME'])\
                                     + linesep + ' '*4
    getset[RUN_READ_CASE_POS] += get_result_file(module) + linesep + ' '*4

    my_study = STUDY_TEMPLATE.format(
        inputs=", ".join(inputs),
        input_names=", ".join(input_names),
        outputs=", ".join(output_names),
        module=module,
        steering_file=steering_file,
        user_fortran=user_fortran,
        run_config_getset=getset[RUN_SET_CONFIG_POS],
        run_read_case_getset=getset[RUN_READ_CASE_POS],
        run_init_getset=getset[RUN_INIT_POS] + getset[RUN_ALLOCATION_POS],
        run_timestep_getset=getset[RUN_TIMESTEP_POS],
        run_finalize_getset=getset[RUN_FINALIZE_POS],
        working_dir=jdc['WORKING_DIRECTORY'],
        res_path=res_folder)

    return my_study


def generate_study_yacs(jdc):
    """
    Creates a yacs file from an eficas dictionary

    @param jdc The eficas data

    @returns A YACS scheme
    """
    try:
        import SALOMERuntime
    except BaseException:
        raise TelemacException(
            "Missing Salome this function "
            "is only available within Salome")

    # Generating python script for the study
    python_script = generate_study_script(jdc)

    # Computing information on inpout/output variables
    input_names, input_types, input_val, output_names, output_types, _ = \
        handle_variables(jdc)

    SALOMERuntime.RuntimeSALOME.setRuntime()
    root = SALOMERuntime.getSALOMERuntime()
    # Building first Bloc
    proc = root.createProc("Function_G")
    proc.setProperty("DefaultStudyID", "1")
    # Adding types
    type_dble = proc.createType("double", "double")
    type_int = proc.createType("int", "int")
    type_string = proc.createType("string", "string")
    type_bool = proc.createType("bool", "bool")
    # Defining container
    cont = proc.createContainer("container0", "Salome")
    cont.clearProperties()

    # Creating function node
    node0 = root.createScriptNode("", "STUDY")
    proc.edAddChild(node0)
    node0.setExecutionMode("remote")
    node0.setScript(python_script)
    node0.setContainer(cont)

    # Adding inputs
    for name, typ, val in zip(input_names, input_types, input_val):
        # Identify yacs type and converting default value
        if typ == b"DOUBLE":
            vartyp = type_dble
            varval = float(val)
        elif typ == b"STRING":
            vartyp = type_string
            varval = val
        elif typ == b"INTEGER":
            vartyp = type_int
            varval = int(val)
        elif typ == b"BOOLEAN":
            vartyp = type_bool
            varval = bool(val)
        else:
            raise TelemacException("Unknow type %s for %s" % (typ, name))
        var = node0.edAddInputPort(name, vartyp)
        var.edInitPy(varval)

    # Adding outputs
    for name, typ in zip(output_names, output_types):
        # Identify yacs type
        if typ == b"DOUBLE":
            vartyp = type_dble
        elif typ == b"STRING":
            vartyp = type_string
        elif typ == b"INTEGER":
            vartyp = type_int
        elif typ == b"BOOLEAN":
            vartyp = type_dble
        else:
            raise TelemacException("Unknow type %s for %s" % (typ, name))
        node0.edAddOutputPort(name, vartyp)

    return proc


def cat_into(out, filenames):
    """
    Merge first line of the first file in filenames
    and second line for all of them

    @param out Merged file
    @param filenames List of files to be merged
    """
    with open(out, 'w') as fobj:
        for filename in filenames:
            fin = open(filename, "r")
            var = fin.readlines()
            # Writing first line of the first file
            if filename == filenames[0]:
                fobj.write(var[0])
            # Last line of each file
            fobj.write(var[-1]+'\n')
            fin.close()


def parametric_study_post(file_to_consider):
    """
    Will gather results from a parametric study (input, output, result_file)

    @param file_to_consider (list) List of type of files to handle
    """
    for _, filetype in enumerate(file_to_consider):
        # Treatment for data files (ascii file containing parameters value)
        if filetype in ['input', 'output']:
            if path.exists(getcwd()+sep+filetype+'.txt'):
                remove(getcwd()+sep+filetype+'.txt')
            files = glob.glob(filetype+'_*')
            files.sort()
            cat_into(filetype+'.txt', files)
            for filename in files:
                remove(filename)
        else:
            # Treatment for results files files (binairy file)
            root, ext = path.splitext(filetype)
            files = glob.glob(root+'_*'+ext)
            files.sort()
            for iindice, filename in enumerate(files):
                rename(filename, root+'_'+str(iindice)+ext)


if __name__ == "__main__":

    with open('test.comm') as jdcfile:
        JDC = jdcfile.read()
    PARAMS = jdc_to_dict(JDC, ['TELEMAC2D', '_F'])

    MY_STUDY = generate_study_script(PARAMS)

    with open('my_study.py', 'w') as FOBJ:
        FOBJ.write(MY_STUDY)

    MY_YACS = generate_study_yacs(PARAMS)
    MY_YACS.saveSchema('my_study.xml')
