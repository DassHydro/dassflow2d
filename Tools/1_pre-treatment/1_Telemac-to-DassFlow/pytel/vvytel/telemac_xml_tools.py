r"""@author TELEMAC-MASCARET Consortium

    @note ... this work is based on a collaborative effort between
                 TELEMAC-MASCARET consortium members
    @brief
"""
# _____                                    __________________________________
# ____/ Specific TELEMAC Toolbox /_________________________________/
#
# _____                             _______________________________________
# ____/ General XML Toolbox /______________________________________/
#
#    Global dictionnaries to avoid having to read these more than once
#    The keys are the full path to the dictionnaries and therefore
#        allows for <root> and <version> to change
from os import path, mkdir, listdir, environ, utime
from utils.files import copy_file
from utils.exceptions import TelemacException
from postel.parser_output import get_latest_output_files

def copy_cas(cas_file, new_cas_file, module, modifs):
    """
    Creates a new cas file from an existing one modifying some values
    The file will have the same metadata as the old one

    @param cas_file (string) Name of the orginal cas
    @param new_cas_file (string) Name of the new cas
    @param module (string) Name of the telemac-mascaret module
    @param modifs (dict) dict where the key is the keyword
                  to modify and value the one to apply
    """
    from execution.telemac_cas import TelemacCas
    from datetime import datetime
    import time
    dico_file = path.join(environ['HOMETEL'],
                          'sources',
                          module,
                          module+'.dico')
    cas = TelemacCas(cas_file, dico_file)

    # update of cas file
    for key, value in modifs.items():
        cas.set(key, value)

    # Writing new cas file
    cas.write(new_cas_file)

    # Transfering modified data of cas_file to the new_file
    # This is to avoid validate_telemac.py to rerun the validation every time
    date = datetime.fromtimestamp(int(path.getmtime(cas_file)))

    mod_time = time.mktime(date.timetuple())
    utime(new_cas_file, (mod_time, mod_time))

    print('    ~> Creating "{}"'.format(path.basename(new_cas_file)))

    return [], [], [], []

###############################################################################
def get_xml_keys(xml, todo):
    """
        Will read the xml's XML keys based on the template todo.
        +: those with None are must have
        +: those without None are optional and reset if there
        Will add extra keys even if it does ont know what to todo.
    """
    xcpt = []  # try all keys for full report
    done = todo.copy()  # shallow copy is here sufficient
    for key in done:
        if key not in xml.keys():
            if done[key] is None:
                xcpt.append({'name': 'get_xml_keys',
                             'msg': 'cannot find the key: ' + key})
        else:
            done[key] = xml.attrib[key]
    if xcpt != []:
        raise TelemacException(xcpt)  # raise full report
    for key in xml.keys():
        if key not in done:
            done[key] = xml.attrib[key]
    return done
###############################################################################
def copy_file_to_valid_dir(cas, valid_dir):
    """
    Copy all the input files from the TelemacCas cas into valid_dir

    @param cas (TelemacCas) Steering file structure
    @param valid_dir (string) The validation directory
    """
    valid_cas_file = path.join(valid_dir, path.basename(cas.file_name))
    # Writing steering file in valid_dir (additionall keyword could have been
    # added)
    cas.write(valid_cas_file)

    cas_dir = path.dirname(cas.file_name)
    # ~~> process sortie files if any
    sortie_files = get_latest_output_files(valid_cas_file)
    # ~~> Copying input files
    for k in cas.in_files:
        submit = cas.in_files[k].split(';')
        val = cas.values[k]
        # Skipping steering file (we already did it)
        if submit[-1] == 'CAS':
            continue
        if val == '':
            continue
        if path.isdir(path.join(cas_dir, val)):
            if not path.exists(path.join(valid_dir, val)):
                mkdir(path.join(valid_dir, val))
            for fle in listdir(path.join(cas_dir, val)):
                if fle == val:
                    continue
                copy_file(path.join(cas_dir, val, fle),
                          path.join(valid_dir, val))
        else:
            copy_file(path.join(cas_dir, val), valid_dir)

    return sortie_files
###############################################################################
def find_targets(dido, src):
    """
    Identify string (in format ref:tag) from the actions.
    For example (vnv_1:T2DRES) will return the result file from action vnv_1

    @param dido (dict) Contains xml actions informations
    @param src (string) The file to find
    """
    layer = []

    if src in dido:
        if src == 'outrefs':
            layer = [[], '', src]
            for k in dido[src]:
                layer[0].append(path.join(dido['path'], k, dido[src][k]))
        else:
            layer = [dido[src], '', src]
    if layer == [] and 'input' in dido:
        for k in dido['input']:
            submit = dido['input'][k].split(';')
            if src in submit[1]:  # filename, filefrom, filetype
                # /!\ Temporary fix because TOMAWAC's IOs names are not yet
                # standard TELEMAC
                if submit[5] == 'SCAL':
                    submit[5] = submit[1]
                # \!/
                layer = [[path.join(dido['safe'], dido['cas'].values[k])],
                         submit[3], submit[5]]
    if layer == [] and 'output' in dido:
        if dido['code'] == 'postel3d':
            for k in dido['output']:
                fle = dido['cas'].values[k]
                if src == path.basename(fle):
                    layer = [[path.join(dido['safe'], fle)], '', 'SELAFIN']
        else:
            for k in dido['output']:
                submit = dido['output'][k].split(';')
                if src in submit[1]:  # filename, filefrom, filetype
                    # /!\ Temporary fix because TOMAWAC's IOs names are not yet
                    # standard TELEMAC
                    if submit[5] == 'SCAL':
                        submit[5] = submit[1]
                    # \!/
                    layer = [[path.join(dido['safe'], dido['cas'].values[k])],
                             submit[3], submit[5]]
    if layer == [] and 'links' in dido:
        for mod in dido['links']:
            if layer == [] and 'in_files' in dido['links'][mod]:
                for k in dido['links'][mod]['input']:
                    submit = dido['links'][mod]['input'][k].split(';')
                    if src in submit[1]:  # filename, filefrom, filetype
                        # /!\ Temporary fix because TOMAWAC's IOs names are not
                        # yet standard TELEMAC
                        if submit[5] == 'SCAL':
                            submit[5] = submit[1]
                        # \!/
                        fle = dido['links'][mod]['cas'].values[k]
                        layer = [[path.join(dido['safe'], fle)],
                                 submit[3], submit[5]]
            if layer == [] and 'output' in dido['links'][mod]:
                if dido['links'][mod]['code'] == 'postel3d':
                    for k in dido['links'][mod]['output']:
                        fle = dido['links'][mod]['cas'].values[k]
                        if src == path.basename(fle):
                            layer = [[path.join(dido['safe'], fle)],
                                     '', 'SELAFIN']
                else:
                    for k in dido['links'][mod]['output']:
                        submit = dido['links'][mod]['output'][k].split(';')
                        if src in submit[1]:  # filename, filefrom, filetype
                            # /!\ Temporary fix because TOMAWAC's IOs names are
                            # not yet standard TELEMAC
                            if submit[5] == 'SCAL':
                                submit[5] = submit[1]
                            # \!/
                            fle = dido['links'][mod]['cas'].values[k]
                            layer = [[path.join(dido['safe'], fle)],
                                     submit[3], submit[5]]

    if layer == [] and 'type' in dido:
        if dido['type'] == src:
            layer = [[dido['target']], '', src]
    return layer
