"""@author Sebastien E. Bourban and Noemie Durand

   @brief All the functions around the compilation
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import re
import shutil
from time import sleep
from os import path, sep, walk, chdir, remove, mkdir, \
    listdir, getcwd
from subprocess import STDOUT, check_output, CalledProcessError
from multiprocessing.sharedctypes import Value, Array
# ~~> dependencies towards other pytel/modules
from utils.files import create_directories, put_file_content, is_newer
from utils.messages import Messages
from utils.progressbar import ProgressBar
from utils.exceptions import TelemacException
from compilation.parser_fortran import get_principal_wrap_names, \
                                       refactor_sources, \
                                       trim_tree, scan_sources
from compilation.scan_tools import get_scan_content, put_scan_content
from config import CFGS

# _____                  ___________________________________________
# ____/ General Toolbox /__________________________________________/
#
LIST_LIBS = {'api':['api', \
    'telemac3d', \
    'telemac2d', \
    'sisyphe', \
    'gaia', \
    'nestor', \
    'tomawac', \
    'artemis', \
    'khione', \
    'waqtel', \
    'utils|gretel', \
    'utils|partel', \
    'utils|bief', \
    'utils|parallel', \
    'utils|hermes', \
    'utils|damocles', \
    'utils|special', \
    'mascaret', \
            ],
             'hermes':[\
    'utils|hermes', \
    'utils|special', \
    ]
            }



def compile_mascaret_dependencies(rebuild, homeres):
    """
    Compile the c files needed by mascaret

    @param rebuild (boolean) If true force to recompile
    @param homeres (dict) Tree dependencies of main program
    """
    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]

    # We need this command to compile c files
    if cfg['cmd_obj_c'] == '':
        raise TelemacException(\
             "\nMissing cmd_obj_c in your configuration file "
             "to compile mascaret C dependencies")

    src_names = []
    src_names.append('Deriv|adstack.c')
    src_names.append('API|cpp|apimascaret.c')

    mes = Messages(size=10)
    # Loop on all of the C files of mascaret
    for src_name in src_names:

        obj_name = src_name.split('|')[-1].replace('.c', '.o')

        src_path = path.join(cfg['root'], 'sources', 'mascaret', \
            *src_name.split('|'))
        obj_path = path.join(cfg['root'], 'builds', cfgname, 'obj', \
            'mascaret', obj_name)

        cmd = cfg['cmd_obj_c'].replace('<srcName>', src_path)
        cmd = cmd.replace('<objName>', obj_path)

        if not path.exists(src_path):
            raise TelemacException( \
                'Could not find the following file' \
                'for compilation: '+path.basename(src_path)+'\n            ' \
                '... so it may have to be removed from the '\
                'following cmdf file: mascaret.cmdf' \
                                  )
        if (is_newer(src_path, obj_path) == 1) and rebuild < 2:
            print('        +> There is no need to compile C object')
        else:
            _, code = mes.run_cmd(cmd, False)
            if code != 0:
                raise TelemacException(\
                        'Could not compile your file adstack')

            print('    - completed: .../sources/mascaret/'+ \
                  src_name.replace('|', '/'))

        homeres['HOMERE_MASCARET']['add'].append((path.dirname(obj_name),
                                                  obj_name, 'mascaret'))

def create_obj_files(oname, oprog, odict, mes, tasks, bypass, homeres, verbose):
    """
    Start the process to compile obj files

    @param oname (string) Name of the file to compile
    @param oprog (string) Name of the main program
    @param odict (dict) Information on file (path, type, libname)
    @param mes (Messages) Structure for execution
    @param tasks (list) Pool of process to run
    @param bypass (boolean) If True bypass errors
    @param homeres (dict) Tree dependencies of the main program
    @param verbose (bool) If true display command

    @returns True if nothign was done (i.e. file up to date)
    """
    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]
    # ~~ Assumes that the source filenames are in lower case ~~~~~~~~
    root, _ = path.splitext(path.basename(oname))

    # ~~ Taggings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    lib_name = odict['libname'].split('.')[0]
    if homeres[oprog]['lib'] != '':
        lib_name = homeres[oprog]['lib']

    # ~~ Directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj_dir = odict['path'].replace(path.join(cfg['root'], 'sources'), \
                                   path.join(cfg['root'], 'builds',
                                             cfgname, 'obj'))
    create_directories(obj_dir)
    chdir(obj_dir)

    # ~~ Removes existing objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if odict['type'][0] == 'M':
        mod_file = path.join(obj_dir, odict['type'][1] + \
                                      cfg['SYSTEM']['sfx_mod'])
        if path.exists(mod_file):
            remove(mod_file)
    obj_file = path.join(obj_dir, root + cfg['SYSTEM']['sfx_obj'])
    if path.exists(obj_file):
        remove(obj_file)

    # ~~ creation of the module:
    cmd = cfg['MODULES'][lib_name]['xobj']
    incs = cfg['MODULES'][lib_name]['incs']
    cmd = cmd.replace('<incs>', incs)
    mods = ''
    if homeres[oprog]['lib'] != '':
        for ones in homeres[oprog]['deps']:
            mod = ones.split('.')[0]
            mods += cfg['MODULES'][lib_name]['mods'] \
                          .replace('<config>', \
                path.join(cfg['MODULES'][mod]['path'], \
                '.'+homeres[oprog]['lib'])) \
                .replace(path.join(cfg['root'], 'sources'), \
                path.join(cfg['root'], 'builds', cfgname, 'obj')) \
                + ' '
    else:
        for ones in homeres[oprog]['deps']:
            mod = ones.split('.')[0]
            mods = mods + \
                cfg['MODULES'][lib_name]['mods'] \
                .replace('<config>', cfg['MODULES'][mod]['path']) \
                .replace(path.join(cfg['root']+sep+'sources'), \
                path.join(cfg['root'], 'builds', cfgname, 'obj')) + ' '

    cmd = cmd.replace('<mods>', mods)
    cmd = cmd.replace('<f95name>', path.join(odict['path'], oname))
    cmd = cmd.replace('<objname>', path.join(odict['path'], obj_file))
    cmd = cmd.replace('<config>', obj_dir).replace('<root>', cfg['root'])

    if verbose:
        print(cmd)
    # ~~> remove ghosts
    out = mes.clean_cmd(tasks)
    _ = mes.start_cmd(tasks, \
        (cmd, bypass, Array('c', b' '*10000), Value('i', 0)), \
        path.join(odict['path'], oname) \
        .replace(path.dirname(cfg['root']), '...'))
    if odict['type'][0] == 'M':
        out.extend(mes.flush_cmd(tasks))
    # ~~> and remove .f from obj_list
    odict['time'] = 1
    return out

def create_lib_files(lname, lmdul, lprog, mprog, mes, tasks,
                     bypass, homeres, verbose):
    """
    Starts the process to compile a library

    @param lname (string) Name of the library
    @param lmdul (string) Name of the module containing the library
    @param lprog (string) Name of the main program associated
    @param mprog (string) Name of the module of the main program
    @param mes (Messages) Structure for execution
    @param tasks (list) Pool of process to run
    @param bypass (boolean) If True bypass errors
    @param homeres (dict) Tree dependencies of the main program
    @param verbose (bool) If true display command
    """
    # ~~ Assumes that all objects are in <config> ~~~~~~~~~~~~~~~~~~~
    # /!\ why path.basename(lname) ?
    # lname takes the successive values of homeres[lprog]['deps']
    #     so it can be passive.solve
    # lmdul is passive
    # lprog is passive.solve
    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]

    # ~~ Taggings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    lib_name = lname.split('.')[0]
    if lname == lprog:
        lib_name = lmdul

    # ~~ May not wish to go the extra mile ~~~~~~~~~~~~~~~~~~~~~~~~~~
    cmd = cfg['MODULES'][lmdul]['xlib']
    if cmd == '':
        return True

    # ~~ Directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj_dir = cfg['MODULES'][lib_name]['path'] \
        .replace(path.join(cfg['root'], 'sources'), \
        path.join(cfg['root'], 'builds', cfgname, 'obj'))
    chdir(obj_dir)
    lib_dir = path.join(cfg['root'], 'builds', cfgname, 'lib')
    create_directories(lib_dir)

    # ~~ Lists all dependent libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dynamic linking
    if 'dyn' in cfg['options']:
        #Path containing all the libraries
        lib_dir = path.join(cfg['root'], 'builds', cfgname, 'lib')
        lib_ext = cfg['SYSTEM']['sfx_lib']

        lib_files = ' -L'+lib_dir+' '
        for lib in homeres[lprog]['deps'][:homeres[lprog]['deps'].index(lname)]:
            l = lib+'4'+mprog
            if not path.exists(path.join(lib_dir, 'lib'+l+lib_ext)):
                raise TelemacException('\nLibrary missing:\n        '+l)
            lib_files += '-l' + l + ' '

    #Static linking
    else:
        lib_files = ''
        for lib in homeres[lprog]['deps'][:homeres[lprog]['deps'].index(lname)]:
            l = path.join(cfg['root'], 'builds', cfgname, 'lib', \
                'lib'+lib+'4'+mprog+cfg['SYSTEM']['sfx_lib'])
            if not path.exists(l):
                raise TelemacException('\nLibrary missing:\n        '+l)
            lib_files = l + ' ' + lib_files

    # ~~ Add external libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if 'libs' in cfg['MODULES'][lmdul]:
        lib_files = lib_files + cfg['MODULES'][lmdul]['libs']

    # lib_file is now created directly within prg[0]'s directory
    # - /!\ hopefuly, the directory exists
    if lmdul == lib_name:
        lib_file = path.join(cfg['root'], 'builds', cfgname, 'lib', \
            'lib'+mprog + cfg['SYSTEM']['sfx_lib'])
    else:
        lib_file = path.join(cfg['root'], 'builds', cfgname, 'lib', \
            'lib'+lib_name+'4'+mprog + cfg['SYSTEM']['sfx_lib'])

    # ~~ Lists all objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj_files = ''
    for pth, obj, lib in homeres[lprog]['add']:
        obj = path.splitext(path.join(pth, obj))[0] \
            .replace(path.join(cfg['root'], 'sources'), \
            path.join(cfg['root'], 'builds', cfgname, 'obj')) \
            .replace(obj_dir, '.') \
            + cfg['SYSTEM']['sfx_obj']
        if lib == lname:
            obj_files = obj_files + (obj+' ')
    for pth, obj, lib in homeres[lprog]['tag']:
        obj = path.splitext(path.join(pth, obj))[0] \
            .replace(path.join(cfg['root'], 'sources'), \
            path.join(cfg['root'], 'builds', cfgname, 'obj')) \
            .replace(obj_dir, '.') \
            +cfg['SYSTEM']['sfx_obj']
        if lib == lname:
            obj_files = obj_files + obj + ' '

    # ~~ is linkage necessary ? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if cfg['COMPILER']['REBUILD'] > 0 and \
        cfg['COMPILER']['REBUILD'] < 2 \
        and path.exists(lib_file):
        remove(lib_file)

    if path.exists(lib_file):
        refresh = False
        for obj in obj_files.split():
            refresh = refresh or (is_newer(obj, lib_file) == 0) \
                or (not path.exists(obj))
        if refresh:
            remove(lib_file)
    if path.exists(lib_file):
        return True

    # ~~ creation of the librairies
    cmd = cmd.replace('<libs>', lib_files)
    if lprog == 'HOMERE_MASCARET':
        cmd = cmd.replace('<objs>', '*.o')
    else:
        cmd = cmd.replace('<objs>', obj_files)
    cmd = cmd.replace('<libname>', lib_file)

    if verbose:
        print(cmd)
    mes.start_cmd(tasks, \
        (cmd, bypass, Array('c', b' '*10000), Value('i', 0)), \
        lib_file.replace(path.dirname(cfg['root']), '...'))
    return False

def create_exe_files(ename, emdul, eprog, mes, bypass, homeres, verbose):
    """
    Compile an executable

    @param ename (string) Name of the program (index in homeres)
    @param emdul (string) Name of the module containing the program
    @param eprog (string) Name of the program to create
    @param mes (Messages) Structure for execution
    @param bypass (boolean) If True bypass errors
    @param homeres (dict) Tree dependencies of the main program
    @param verbose (bool) If true display command

    @returns True if nothign was done (i.e. file up to date)
    """

    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]
    # ~~ Directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    lib_dir = path.join(cfg['root'], 'builds', cfgname, 'lib')
    obj_dir = cfg['MODULES'][emdul]['path'] \
        .replace(path.join(cfg['root'], 'sources'), \
        path.join(cfg['root'], 'builds', cfgname, 'obj'))
    # this is because of the aggregation of the local objects
    # within the executable project
    chdir(obj_dir)
    exe_dir = path.join(cfg['root'], 'builds', cfgname, 'bin')
    create_directories(exe_dir)

    # ~~ Command line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cmd = cfg['MODULES'][emdul]['xexe']

    # ~~ Executables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    exe_file = path.join(exe_dir, eprog + cfg['SYSTEM']['sfx_exe'])
    if '<exename>' in cfg['MODULES'][emdul]['xexe']:
        cmd = cmd.replace('<exename>', exe_file)\
            .replace('<config>', lib_dir)\
            .replace('<root>', cfg['root'])
    obj_cmd = path.join(obj_dir, eprog + '.cmdo')
    exe_cmd = path.join(lib_dir, eprog + '.cmdx')

    # ~~ Lists all system libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    lib_files = ''
    # /!\ [1:] to create the exe from local objs.
    if 'dyn' in cfg['options']:
        #Path containing all the libraries
        lib_dir = path.join(cfg['root'], 'builds', cfgname, 'lib')
        lib_ext = cfg['SYSTEM']['sfx_lib']

        lib_files = ' -L'+lib_dir+' '
        for lib in homeres[ename]['deps'][:-1]:
            l = lib+'4'+eprog
            if not path.exists(path.join(lib_dir, 'lib'+l+lib_ext)):
                raise TelemacException('\nLibrary missing:\n        '+l)
            lib_files += '-l' + l + ' '

    #Static linking
    else:
        lib_files = ''
        for lib in homeres[ename]['deps'][:-1]:
            l = path.join(cfg['root'], 'builds', cfgname, 'lib', \
                'lib'+lib+'4'+eprog+cfg['SYSTEM']['sfx_lib'])
            if not path.exists(l):
                raise TelemacException('\nLibrary missing:\n        '+l)
            lib_files = l + ' ' + lib_files

    for lib in cfg['ADDONES']:
        for fli in cfg['ADDONES'][lib]:
            ones = lib+'.'+path.splitext(fli)[0]
            if not ones in homeres:
                raise TelemacException(\
                    '\nYou may have forgotten to rescan for: ' \
                    +ones+' not in homeres (' \
                    +repr(homeres.keys()))
            for one in homeres[ones]['deps'][:len(homeres[ones]['deps'])-1]:
                l = path.join(cfg['root'], 'builds', cfgname, 'lib', \
                    one+'4'+ones+cfg['SYSTEM']['sfx_lib'])
                if not path.exists(l):
                    raise TelemacException(\
                        '\nAdded library missing:\n        '+l)
                lib_files = l + ' ' + lib_files
            l = path.join(cfg['root'], 'builds', cfgname, 'lib', \
                ones+cfg['SYSTEM']['sfx_lib'])
            if not path.exists(l):
                raise TelemacException('\nLibrary missing:\n        '+l)
            lib_files = l + ' ' + lib_files
    lib = homeres[ename]['deps'][len(homeres[ename]['deps'])-1]
    if lib == emdul:
        lib_file = path.join(cfg['root'], 'builds', cfgname, 'lib', \
            'lib'+eprog+cfg['SYSTEM']['sfx_lib'])
    else:
        lib_file = path.join(cfg['root'], 'builds', cfgname, 'lib', \
            'lib'+lib+'4'+eprog+cfg['SYSTEM']['sfx_lib'])
    if not path.exists(lib_file):
        raise TelemacException('\nLibrary missing:\n        '+lib_file)

    # ~~ Add external libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if 'libs' in cfg['MODULES'][emdul]:
        lib_files = lib_files + cfg['MODULES'][emdul]['libs']

    # ~~ Lists local objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj_files = ''
    for pth, obj, lib in homeres[ename]['add']:
        root, _ = path.splitext(path.basename(path.join(pth, obj)))
        obj_name = path.basename(obj).lower()+cfg['SYSTEM']['sfx_obj']
        if lib == emdul and obj_name  not in obj_files.split():
            obj_file = root.lower()+cfg['SYSTEM']['sfx_obj']
            if not path.exists(obj_file):
                raise TelemacException('Object missing:\n        '+obj_file)
            obj_files = obj_files + obj_file + ' '

    for pth, obj, lib in homeres[ename]['tag']:
        root, _ = path.splitext(path.basename(path.join(pth, obj)))
        obj_name = path.basename(obj).lower()+cfg['SYSTEM']['sfx_obj']
        if lib == emdul and obj_name not in obj_files.split():
            obj_file = root.lower()+cfg['SYSTEM']['sfx_obj']
            if not path.exists(obj_file):
                raise TelemacException('Object missing:\n        '+obj_file)
            obj_files = obj_files + obj_file +' '

    # ~~ is executable necessary ? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if path.exists(exe_file):
        if cfg['COMPILER']['REBUILD'] > 0 and cfg['COMPILER']['REBUILD'] < 3:
            remove(exe_file)
    if path.exists(exe_file):
        if cfg['COMPILER']['REBUILD'] > 2 or cfg['COMPILER']['REBUILD'] == 0:
            refresh = False
            for obj_file in obj_files.split():
                refresh = refresh or (is_newer(obj_file, exe_file) == 0)
            for l in lib_files.split():
                # Only checks the telemac libraries
                if l.find(path.join(cfg['root'], 'builds',
                                    cfgname, 'lib')) != -1:
                    refresh = refresh or (is_newer(l, exe_file) == 0)
            if refresh:
                remove(exe_file)
    if path.exists(exe_file):
        return True

    # ~~ creation of the exe (according to makefile.wnt + systel.ini):
    cmd = cmd.replace('<libs>', lib_files)
    if ename == 'HOMERE_MASCARET':
        cmd = cmd.replace('<objs>', '*.o')
    else:
        cmd = cmd.replace('<objs>', obj_files)

    xocmd = cfg['MODULES'][emdul]['xobj']
    xocmd = xocmd.replace('<incs>', cfg['MODULES'][emdul]['incs'])
    mods = ''
    for mod in homeres[ename]['deps']:
        mods = mods + cfg['MODULES'][emdul]['mods'] \
            .replace('<config>', cfg['MODULES'][mod]['path']) \
            .replace(path.join(cfg['root'], 'sources'), \
            path.join(cfg['root'], 'builds', cfgname, 'obj')) \
            + ' '
    xocmd = xocmd.replace('<mods>', mods)
    # <> ... still to be replaced
    xocmd = xocmd.replace('<config>', lib_dir) \
        .replace('<root>', cfg['root'])

    lib_files = lib_file + ' ' + lib_files
    xecmd = cfg['MODULES'][emdul]['xexe']
    xecmd = xecmd.replace('<libs>', lib_files)
    # Special keyword for nag with ', ' separating the libraries
    xecmd = xecmd.replace('<libsnag>', lib_files.replace(' ', ','))
    # <exename> and <objs> ... still to be replaced
    xecmd = xecmd.replace('<config>', lib_dir).replace('<root>', cfg['root'])

    if verbose:
        print(cmd)
    tail, code = mes.run_cmd(cmd, bypass)
    if tail != b'':
        if path.exists(obj_cmd):
            remove(obj_cmd)
        if path.exists(exe_cmd):
            remove(exe_cmd)
    if code != 0:
        if path.exists(obj_cmd):
            remove(obj_cmd)
        if path.exists(exe_cmd):
            remove(exe_cmd)
        raise TelemacException(\
                'Could not link your executable. ' \
                'Please verify your external library installation' \
                ' or the python script itself.')
    print('    - created ' + exe_file.replace(path.dirname(cfg['root']), '...'))

    # ~~> Make the keys portable (no full path)
    for k in cfg['TRACE']:
        xocmd = xocmd.replace(cfg['TRACE'][k], '['+k+']')
        xecmd = xecmd.replace(cfg['TRACE'][k], '['+k+']')
    put_file_content(obj_cmd, [xocmd])
    put_file_content(exe_cmd, [xecmd])

    return False

def create_pyd_files(yname, yfile, ymdul, yprog, mes, bypass, homeres, verbose):
    """
    Generates pyd files (alternate method to compile_api)

    @param yname Name of the program (homeres index)
    @param yfile Name of the file
    @param ymdul Name of the module
    @param yprog Name of the program
    @param mes (Messages) Structure for execution
    @param bypass (boolean) If True bypass errors
    @param homeres (dict) Tree dependencies of the main program
    @param verbose (bool) If true display command
    """

    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]
    # ~~ Directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    lib_dir = path.join(cfg['root'], 'builds', cfgname, 'lib')
    pyd_dir = path.join(cfg['root'], 'builds', cfgname, 'bin')
    create_directories(pyd_dir)
    chdir(pyd_dir)

    # ~~ Command line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cmd = cfg['MODULES'][ymdul]['xpyd']

    # ~~ Python module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if 'sfx_pyd' not in cfg['SYSTEM']:
        raise TelemacException(\
            'Suffix extension missing from your ' \
            'configuration file (just add sfx_pyd:.pyd)')
    pyd_file = path.join(pyd_dir, yprog + cfg['SYSTEM']['sfx_pyd'])
    if '<pydname>' in cfg['MODULES'][ymdul]['xpyd']:
        cmd = cmd.replace('<pydname>', yprog) \
            .replace('<config>', lib_dir) \
            .replace('<root>', cfg['root'])
    # ~~~
    f95_file = path.join(cfg['MODULES'][ymdul]['path'], yfile)
    if '<>' in cfg['MODULES'][ymdul]['xpyd']:
        cmd = cmd.replace('<f95name>', f95_file)

    # ~~ Lists all system libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    lib_files = ''
    # /!\ [1:] to create the exe from local objs.
    for lib in homeres[yname]['deps'][:len(homeres[yname]['deps'])-1]:
        l = path.join(cfg['root'], 'builds', cfgname, 'lib', \
            lib+'4'+yprog+cfg['SYSTEM']['sfx_lib'])
        if not path.exists(l):
            raise TelemacException('Library missing:\n        '+l)
        lib_files = l + ' ' + lib_files
    for lib in cfg['ADDONES']:
        for fli in cfg['ADDONES'][lib]:
            ones = lib+'.'+path.splitext(fli)[0]
            if not ones in homeres:
                raise TelemacException(\
                    'You may have forgotten to rescan for: ' + \
                    ones+' not in homeres (' + \
                    repr(homeres.keys()))
            for one in homeres[ones]['deps'][:len(homeres[ones]['deps'])-1]:
                l = path.join(cfg['root'], 'builds', cfgname, 'lib', \
                    one+'4'+ones+cfg['SYSTEM']['sfx_lib'])
                if not path.exists(l):
                    raise TelemacException(\
                        'Added library missing:\n        '+l)
                lib_files = l + ' ' + lib_files
            l = path.join(cfg['root'], 'builds', cfgname, 'lib', \
                ones+cfg['SYSTEM']['sfx_lib'])
            if not path.exists(l):
                raise TelemacException('Library missing:\n        '+l)
            lib_files = l + ' ' + lib_files
    lib = homeres[yname]['deps'][len(homeres[yname]['deps'])-1].split('.')[0]
    lib_file = path.join(cfg['root'], 'builds', cfgname, 'lib', \
        yprog+cfg['SYSTEM']['sfx_lib'])
    if not path.exists(lib_file):
        raise TelemacException('Library missing:\n        '+lib_file)

    # ~~ Add external libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if 'libs' in cfg['MODULES'][ymdul]:
        lib_files = lib_files + cfg['MODULES'][ymdul]['libs']

    # ~~ Lists local objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj_files = ''
    for obj, lib in homeres[yname]['add']:
        root, _ = path.splitext(path.basename(obj))
        if lib == ymdul and\
            obj.lower()+cfg['SYSTEM']['sfx_obj'] not in obj_files.split():
            obj_file = root.lower()+cfg['SYSTEM']['sfx_obj']
            if not path.exists(obj_file):
                raise TelemacException('Object missing:\n        '+obj_file)
            obj_files = obj_files + obj_file + ' '
    #if obj_files.strip() == '' and path.exists(exe_file): return True
    for obj, lib in homeres[yname]['tag']:
        obj_name = (path.basename(obj)).lower()+cfg['SYSTEM']['sfx_obj']
        if lib == ymdul and obj_name not in obj_files.split():
            if not path.exists(obj_name):
                raise TelemacException('Object missing:\n        '+obj_name)
            obj_files = obj_files + obj_name +' '
    # ~~ Lists local files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # ~~ is executable necessary ? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if path.exists(pyd_file):
        if cfg['COMPILER']['REBUILD'] > 0 and cfg['COMPILER']['REBUILD'] < 3:
            remove(pyd_file)
    if path.exists(pyd_file):
        if cfg['COMPILER']['REBUILD'] > 2 or cfg['COMPILER']['REBUILD'] == 0:
            refresh = False
            for obj_file in obj_files.split():
                refresh = refresh or (is_newer(obj_file, pyd_file) == 0)
            for l in lib_files.split():
                # Only checks the telemac libraries
                if l.find(path.join(cfg['root'], 'builds',
                                    cfgname, 'lib')) != -1:
                    refresh = refresh or (is_newer(l, pyd_file) == 0)
            if refresh:
                remove(pyd_file)
    if path.exists(pyd_file):
        return True

    # ~~ creation of the exe
    cmd = cmd.replace('<libs>', lib_files.strip())
    cmd = cmd.replace('<objs>', obj_files.strip())

    cmd = cmd.replace('<incs>', cfg['MODULES'][ymdul]['incs'].strip())
    mods = ''
    for mod in homeres[yname]['deps'][:-1]:
        mods = cfg['MODULES'][ymdul]['mods'] \
            .replace('<config>', \
            cfg['MODULES'][mod.split('.')[0]]['path'] \
            .replace(path.join(cfg['root'], 'sources'), \
            path.join(cfg['root'], 'builds', cfgname, 'obj'))) \
            + ' ' + mods
    cmd = cmd.replace('<mods>', mods.strip())

    lib_files = lib_file + ' ' + lib_files

    if verbose:
        print(cmd)
    _, code = mes.run_cmd(cmd, False)
    if code != 0:
        raise TelemacException(\
                'Could not link your executable. ' \
                'Please verify your external library installation ' \
                'or the python script itself.')
    print('    - created ' + pyd_file.replace(path.dirname(cfg['root']), '...'))

    return False

def mycopy(src, dst):
    """
    Custom copy that will remove the destination first if it is there
    @param src The file to copy
    @param dst The destiantion where to copy the file
    """
    if path.exists(dst):
        remove(dst)
    shutil.copy(src, dst)

def get_api_incs_flags():
    """
    Retuns the string for incs_flags for api

    @returns the string
    """
    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]
    api_dir = path.join(cfg['root'], 'builds', cfgname, 'wrap_api')

    incs_flags = cfg.get('incs_all', '') \
        .replace('<root>', cfg['root']) \
        .replace('\n', ' ')
    incs_flags += ' -I'+api_dir+sep+'include'

    return incs_flags

def get_api_ld_flags(api_name, static):
    """
    Retuns the string for ld_flags for api

    @param api_name Name of api (hermes ot api)
    @param static If true libraries are considered static

    @returns the string
    """
    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]
    api_dir = path.join(cfg['root'], 'builds', cfgname, 'wrap_api')
    lib_dir = path.join(cfg['root'], 'builds', cfgname, 'lib')
    lib_api_dir = path.join(api_dir, 'lib')

    ld_flags = ''
    # Adding library path
    if not static:
        ld_flags += ' -L'+lib_api_dir + ' '
    else:
        ld_flags += ' '

    lib_ext = cfg['sfx_lib']
    # Adding list of libraries
    for lib_name in LIST_LIBS[api_name]:
        if lib_name == 'mascaret':
            if not path.exists(path.join(lib_dir, 'libmascaret'+lib_ext)):
                continue

        if lib_name in ['mascaret', 'api']:
            # Not adding 4api
            lib = lib_name.split('|')[-1]
        else:
            lib = lib_name.split('|')[-1]+'4api'

        if static:
            ld_flags += lib_dir+sep+'lib'+lib+lib_ext+" "
        else:
            ld_flags += "-l"+lib+" "

    ld_flags += cfg['libs_all'] \
        .replace('<root>', cfg['root']) \
        .replace('\n', ' ')

    return ld_flags

def compile_princi_lib(princi_file, incs_flags, ld_flags):
    """
    Compiling user fortran as a library

    @param princi_file (string) Path of user fortran
    @param incs_flags Include flags for compilation
    @param ld_flags Linking flags for compilation
    """
    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]
    if not path.exists(princi_file):
        raise TelemacException(\
            'Could not find your FORTRAN: '+princi_file)

    user_fortran = []
    # in case of a folder getting list of files
    if path.isdir(princi_file):
        list_files = listdir(princi_file)
        for fle in list_files:
            if re.match("^m[0-9]+.*", fle) and \
                fle.lower().endswith((".f", ".f90")):
                user_fortran.append(path.join(princi_file, fle))
        # Adding the other files
        for fle in list_files:
            if fle not in user_fortran and \
                fle.lower().endswith((".f", ".f90")):
                user_fortran.append(path.join(princi_file, fle))
    else:
        user_fortran = [princi_file]
    # Building linking commands
    command = cfg['cmd_lib'] \
        .replace('<libname>', "libuser_fortran" + cfg['sfx_lib']) \
        .replace('<objs>', ' '.join(user_fortran)) \
        .replace('<libs>', '')
    command += ' ' + incs_flags + ' ' + ld_flags

    mes = Messages(size=10)
    tail, code = mes.run_cmd(command, False)
    if code != 0:
        raise TelemacException(\
            'Could not compile your FORTRAN \
            (runcode='+str(code)+').\n        '+tail)

def copy_src_api(api_dir, src_list, src_dir):
    """
    Copying sources in src_list from src_dir in api folder renaming them in .f90
    And build a string containing list of names

    @param api_dir Path to the api directory
    @param src_list List of files to copy
    @param src_dir Path of the sources

    @returns String containing list of new names for f2py
    """

    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]
    lib_dir = path.join(cfg['root'], 'builds', cfgname, 'lib')

    source = ''
    for src in src_list:
        root, _ = path.splitext(src)
        # Copying source in wrap_api folder and
        # changin extension into .f90
        mycopy(path.join(src_dir, src), \
            path.join(api_dir, 'src', root+'.f90'))
        # Building list of sources
        source += path.join(api_dir, 'src', root+'.f90') + ' '

    # Copying libraries
    dyn_ext = cfg['sfx_lib']
    for lib in LIST_LIBS['api']:
        lib_name = lib.split('|')[-1]
        lib_name_tel = 'lib'+lib_name+'4api'+dyn_ext
        # Mascaret is not named the same
        if lib == 'mascaret':
            lib_name_tel = 'lib'+lib_name+dyn_ext
            if not path.exists(path.join(lib_dir, lib_name_tel)):
                continue
        if lib == 'api':
            lib_name_tel = 'libapi'+dyn_ext
        mycopy(path.join(lib_dir, lib_name_tel), \
            path.join(api_dir, 'lib', lib_name_tel))

    return source

def generate_api():
    """
    Builds the structure for an Python API

    @returns source list for api and for hermes
    """
    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]

    api_dir = path.join(cfg['root'], 'builds', cfgname, 'wrap_api')
    obj_dir = path.join(cfg['root'], 'builds', cfgname, 'obj')
    if not path.exists(api_dir):
        mkdir(api_dir)
        mkdir(api_dir+sep+'lib')
        mkdir(api_dir+sep+'src')
        mkdir(api_dir+sep+'include')

    # Copying Modules
    for root, _, files in walk(obj_dir):
        for ffile in files:
            if ffile.endswith("mod"):
                mycopy(path.join(root, ffile), \
                    path.join(api_dir, 'include', ffile))

    # Copying sources for t3d, t2d, art and sis
    src_list = []
    src_list.append("api_handle_var_t3d.f")
    src_list.append("api_handle_var_t2d.f")
    src_list.append("api_handle_var_sis.f")
    src_list.append("api_handle_var_art.f")
    src_list.append("api_handle_var_wac.f")
    src_list.append("api_handle_error.f")
    src_list.append("api_interface.f")
    src_dir = path.join(cfg['root'], 'sources', 'api')
    source_api = copy_src_api(api_dir, src_list, src_dir)
    # Copying sources for hermes
    src_list = []
    src_list.append("close_bnd.f")
    src_list.append("close_mesh.f")
    src_list.append("get_bnd_connectivity.f")
    src_list.append("get_bnd_ipobo.f")
    src_list.append("get_bnd_nelem.f")
    src_list.append("get_bnd_npoin.f")
    src_list.append("get_bnd_numbering.f")
    src_list.append("get_bnd_value.f")
    src_list.append("get_bnd_color.f")
    src_list.append("get_data_ntimestep.f")
    src_list.append("get_data_nvar.f")
    src_list.append("get_data_time.f")
    src_list.append("get_data_value.f")
    src_list.append("get_data_var_list2.f")
    src_list.append("get_mesh_connectivity.f")
    src_list.append("get_mesh_coord.f")
    src_list.append("get_mesh_date.f")
    src_list.append("get_mesh_dimension.f")
    src_list.append("get_mesh_l2g_numbering.f")
    src_list.append("get_mesh_nelem.f")
    src_list.append("get_mesh_nplan.f")
    src_list.append("get_mesh_npoin.f")
    src_list.append("get_mesh_npoin_per_element.f")
    src_list.append("get_mesh_nptir.f")
    src_list.append("get_mesh_title.f")
    src_list.append("get_mesh_orig.f")
    src_list.append("open_bnd.f")
    src_list.append("open_mesh.f")
    src_list.append("set_bnd.f")
    src_list.append("set_mesh.f")
    src_list.append("set_header.f")
    src_list.append("add_data.f")
    src_list.append("transfer_group_info.f")
    src_list.append("get_error_message.f")
    src_dir = path.join(cfg['root'], 'sources', 'utils', 'hermes')
    source_hermes = copy_src_api(api_dir, src_list, src_dir)

    return source_api, source_hermes

def compile_api_f2py(name, api_dir, source_list, skip_source, \
    ld_flags, f2py_name, fcompiler, compiler, silent, f2py_opt=''):
    """
    Running f2py to generate Python wrapper

    @param name Name of the wrapper
    @param api_dir Path to the api folder
    @param source_list List of source for the api
    @param skip_source List of function to skip
    @param ld_flags Linking flags
    @param f2py_name Name of the f2py executable (f2py by default)
    @param fcompiler Name of the fortran compiler
    @param compiler Name of the c compiler
    @param silent If True f2py is run in silent mode and commad are not
    displayed
    @param f2py_opt Addtional options passed to f2py (option --opt)
    """

    # Generating Py wrapper using f2py
    pyf_file = path.join(api_dir, 'lib', name+'.pyf')
    if path.exists(pyf_file):
        remove(pyf_file)
    if skip_source != '':
        skip_source = 'skip: ' + skip_source + ' :'
    quiet = '--quiet' if silent else ''
    # First step of call to f2py
    cmd = '{f2py_name} {quiet} -h {pyf_file} -m _{name} '\
          '{source_list} {skip_source}' \
        .format(f2py_name=f2py_name,
                quiet=quiet,
                pyf_file=pyf_file,
                name=name,
                source_list=source_list,
                skip_source=skip_source)
    if not silent:
        print(cmd)
    try:
        output = check_output(cmd, shell=True, stderr=STDOUT)
    except CalledProcessError as execpt:
        raise TelemacException(\
                'Error during first part of f2py for {} {} \n'
                '{}'.format(name, execpt.returncode,
                            execpt.output.decode('utf-8')))
    if not silent:
        print(output.decode('utf-8'))
    print("    ~> First part of f2py for {} passed".format(name))

    pwd = getcwd()
    chdir(path.join(api_dir, 'lib'))
    if compiler != '':
        compile_cmd = "--compiler="+compiler
    else:
        compile_cmd = ''

    f2py_build_dir = 'tmp_f2py'

    # Second step of call to f2py
    cmd = '{f2py_name} {quiet} -c {pyf_file} --fcompiler={fcompiler} '\
          '{compile_cmd} --opt="{f2py_opt}" -I{include} {ld_flags} ' \
          '--build-dir {build_dir}' \
        .format(f2py_name=f2py_name,
                pyf_file=pyf_file,
                quiet=quiet,
                fcompiler=fcompiler,
                compile_cmd=compile_cmd,
                include=path.join(api_dir, 'include'),
                ld_flags=ld_flags,
                f2py_opt=f2py_opt,
                build_dir=f2py_build_dir)
    if not silent:
        print(cmd)
    try:
        output = check_output(cmd, shell=True, stderr=STDOUT)
    except CalledProcessError as execpt:
        raise TelemacException(\
                'Error during second part of f2py for {} {} \n'
                '{}'.format(name, execpt.returncode,
                            execpt.output.decode('utf-8')))
    if not silent:
        print(output.decode('utf-8'))
    # Removing build directory
    if path.exists(f2py_build_dir):
        shutil.rmtree(f2py_build_dir)
    print("    ~> Second part of f2py of %s passed"%name)
    chdir(pwd)

def compile_api_files(silent, static=False, hermes_only=False):
    """
    Compiling the APIs for Telemac-Mascaret

    @param silent (boolean) If True does display commands
    @param static If true libraries are considered static
    @param hermes_only If true only the hermes api will be compiled
    """
    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]
    print('\nBuilding the Python API \n'+'~'*72+'\n')

    source_api, source_hermes = generate_api()
    print("    ~> Wrap_api built")

    # Excluding subroutine from api_handle_var_xxx for each module
    skip_source = ''
    for short in ['t2d', 'sis', 't3d', 'art', 'wac']:
        skip_source += 'get_boolean_'+short+'_d get_double_'+short+'_d '
        skip_source += 'get_integer_'+short+'_d get_string_'+short+'_d '
        skip_source += 'get_var_size_'+short+'_d set_boolean_'+short+\
                       '_d set_double_'+short+'_d '
        skip_source += 'set_integer_'+short+'_d set_string_'+short+'_d '
        skip_source += 'get_double_array_'+short+\
                       '_d get_integer_array_'+short+'_d '
        skip_source += 'set_double_array_'+short+\
                       '_d set_integer_array_'+short+'_d '

    api_dir = path.join(cfg['root'], 'builds', cfgname, 'wrap_api')
    compiler = cfg.get('pyd_compiler', '')
    if 'pyd_fcompiler' not in cfg:
        raise TelemacException(\
            "Missing keyword pyd_fcompiler in configuration file")
    fcompiler = cfg['pyd_fcompiler']
    f2py_name = cfg.get('f2py_name', 'f2py')
    f2py_opt = cfg.get("f2py_opt", '')

    if not hermes_only:
        print("    ~> Compiling Modules api")
        ld_flags = get_api_ld_flags('api', static)
        compile_api_f2py('api', api_dir, source_api, skip_source, ld_flags,
                         f2py_name, fcompiler, compiler, silent,
                         f2py_opt=f2py_opt)

    print("    ~> Compiling hermes api")
    ld_flags = get_api_ld_flags('hermes', static)
    compile_api_f2py('hermes', api_dir, source_hermes, '',
                     ld_flags, f2py_name, fcompiler, compiler, silent)

def update_cmdf(bypass, cleanup, verbose):
    """
    Update the cmdf file by doing a scan of the sources

    @param bypass If True Exception are bypassed
    @param cleanup If True cmdf are removed and rewritten
    @param verbose If True will print the scanned file
    """
    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]
    # ~~ Scans all source files to build a relation database ~~~~~~~~~~
    # TODO: parallelistaion of the scan_sources
    fic, _, _, _, _, top, _, whocallswho = \
            scan_sources(cfgname, cfg, bypass, verbose)

    # ~~ Builds the Call Tree for each tree top ~~~~~~~~~~~~~~~~~~~~~~~
    homeres = {}
    print('\nUpdating your cmdf file for compilation without scan\n'
          + '~' * 72 + '\n')
    # ~~> top
    # has the shape of a dictionary file, where the entry key is the
    # name of a main program (e.g. SPLITSEL, HOMERE_PARTEL...)
    # and where the value are single cell array with the name
    # of the lib
    # the main program depends upon (e.g. ['splitsel'], ['partel']...)
    for item in top:
        # TODO: check whether you can expect more than one
        # element in top[item]
        for ones in top[item]:
            mod = ones.split('.')[0]
            # ~~> filtering unwanted modules
            if mod not in cfg['COMPILER']['MODULES']:
                continue
            # ~~> for each of those, you are expected to write a cmdf
            #    file, which will be used to compile that particular
            #    top of the tree and its dependencies.
            print('        +> ' + item + ' (' + ones + ')')

            # ~~ Builds the Call Tree for each main program ~~~~~~~~~~~
            rebuild = cfg['COMPILER']['REBUILD']
            maksystel = {'list': [], 'deps': []}
            maksystel['deps'] = trim_tree(item, ones, whocallswho,
                                          rebuild, maksystel)
            homeres[ones] = maksystel
            # ~~ Prepare the cmdf file to avoid future scans ~~~~~~~~~~
            for_dir = whocallswho[ones]['path']
            if mod in cfg['ADDONES']:
                for_dir = whocallswho[ones]['path']
                for_cmd = path.join(for_dir, ones + '.cmdf')
            elif 'homere' in item.lower() or 'systeme' in item.lower():
                for_cmd = path.join(for_dir,
                                    '_'.join((item.lower())
                                             .split('_')[1:]) + '.cmdf')
            else:
                for_cmd = path.join(for_dir, item.lower() + '.cmdf')
            # /!\ 'name' defines the name of the executable
            exe_name = path.splitext(path.basename(for_cmd))[0]
            file_list = {'general': {'path': whocallswho[ones]['path'],
                                     'name': exe_name,
                                     'module': mod,
                                     'liborder': maksystel['deps']}}
            for obj, lib in homeres[ones]['list']:
                obj_name = path.splitext(path.basename(
                    obj.replace('|', sep)))[0]
                fic = whocallswho[lib].get(obj_name.upper(), None)
                if fic is None:
                    raise TelemacException(\
                                 '+> missmatch between Fortran name '
                                 'and file name for: '
                                 + obj_name.upper())
                if lib not in file_list:
                    file_list.update({lib: {'path': fic['path'],
                                            'files': []}})
                file_list[lib]['files'].append(fic['file'])
            if not path.exists(for_cmd) or rebuild == 2 or cleanup:
                put_scan_content(for_cmd, cfg['root'], file_list)
            else:
                fixe_list = get_scan_content(for_cmd, cfg['root'], bypass)
                # ~~> check the update for new libraries
                if file_list['general']['liborder'] != \
                        fixe_list['general']['liborder']:
                    fixe_list['general']['liborder'] = \
                        file_list['general']['liborder']
                    print(' ' * 9 + 'The number of elements linked '
                                    'together has changed: '
                          + ' | '
                          .join(file_list['general']['liborder']))
                    fixes = fixe_list.keys()
                    for lib in fixes:
                        if lib == 'general':
                            continue
                        if lib not in file_list:
                            del fixe_list[lib]
                    for lib in file_list:
                        if lib == 'general':
                            continue
                        if lib not in fixe_list:
                            fixe_list.update({lib: {
                                'path': file_list[lib]['path'],
                                'files': file_list[lib]['files']}})
                # ~~> add new files
                mes = ''
                for lib in file_list:
                    if lib == 'general':
                        continue
                    if lib in fixe_list:
                        for fic in file_list[lib]['files']:
                            if fic not in fixe_list[lib]['files']:
                                mes += '\n            ~ ' + lib + \
                                       ' | ' + fic
                                fixe_list[lib]['files'].append(fic)
                if mes != '':
                    print(' ' * 9 + 'The following have been added to'
                                    ' the CMDF file: ' + for_cmd + mes)
                # ~~> remove inexistant files
                mes = ''
                for lib in fixe_list:
                    if lib == 'general':
                        continue
                    fixes = fixe_list[lib]['files']
                    for fix in fixes:
                        if not path.exists(
                                path.join(fixe_list[lib]['path'],
                                          fix.replace('|', sep))):
                            mes += '\n            ~ ' + lib + \
                                   ' | ' + fix
                            index = fixe_list[lib]['files'].index(fix)
                            del fixe_list[lib]['files'][index]
                if mes != '':
                    print(' ' * 9 + 'The following will be removed '
                          'from the CMDF file: '
                          + for_cmd + mes)
                # ~~> put content as CMDF file
                put_scan_content(for_cmd, cfg['root'], fixe_list)


def compile_cmdf(ncsize, modules, verbose):
    """
    Compile files for each cmdf

    @param ncsize (integer) Number of processor for parallel compilation
    @param modules (list) If modules is not empty compiling only
                          the modules given in the list
    @param verbose (boolean) If True display compilation command
    """
    # TODO: remove bypass ?
    bypass = False
    cfgname = CFGS.cfgname
    cfg = CFGS.configs[cfgname]
    cmdf_files = {}
    homeres = {}
    found = False
    rebuild = cfg['COMPILER']['REBUILD']
    # Defining list of modules to compile
    if modules != []:
        list_modules = modules
    else:
        list_modules = cfg['COMPILER']['MODULES']
    for mod in sorted(list_modules):
        cmdf_files.update({mod: {}})
        if mod in cfg['MODULES']:
            found = found or (cfg['MODULES'][mod]['cmdfs'] != [])
            # make sure the key cmdfs exists
            for cmdf_file in cfg['MODULES'][mod]['cmdfs']:
                cmdf = get_scan_content(
                    path.join(cfg['MODULES'][mod]['path'], cmdf_file),
                    cfg['root'], bypass)
                cmdf_files[mod].update({cmdf['general']['name']: cmdf})
        # ~~ Look whether .o older than .f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for item in sorted(cmdf_files[mod]):
            print('\n\nCompiling from the tree top ' + item
                  + ' plus dependents\n' + '~' * 72 + '\n')
            maksystel = {'add': [], 'tag': [],
                         'deps': cmdf_files[mod][item]['general']['liborder'],
                         'lib': ''}
            homeres[item] = maksystel
            for lib in maksystel['deps']:
                for_dir = cmdf_files[mod][item][lib]['path']
                if mod in cfg['ADDONES']:
                    homeres[item]['lib'] = mod
                    for_dir = cmdf_files[mod][item][lib]['path']\
                        + sep + '.' + mod
                    create_directories(for_dir)
                for fli in cmdf_files[mod][item][lib]['files']:
                    # In case the file is in a subfolder of the module
                    # replace the | that defines the separator by the os
                    # separator
                    fle = fli.replace('|', sep)
                    # /!\ original file
                    src_name = cmdf_files[mod][item][lib]['path'] + sep + fle
                    obj_dir = for_dir.replace(path.join(cfg['root'], 'sources'),
                                              path.join(cfg['root'], 'builds',
                                                        cfgname, 'obj'))
                    create_directories(obj_dir)
                    obj_name = obj_dir + sep
                    obj_name += path.splitext(path.basename(fle))[0]
                    obj_name = obj_name + cfg['SYSTEM']['sfx_obj']
                    if not path.exists(src_name):
                        raise TelemacException(\
                                '+> Could not find the following file '
                                'for compilation: '
                                + path.basename(src_name)
                                + '\n            ... so it may have to be '
                                'removed from the following cmdf file ')
                    if (is_newer(src_name, obj_name) == 1) \
                         and rebuild < 2:
                        homeres[item]['tag'].append((for_dir, fle, lib))
                    else:
                        if mod in cfg['ADDONES']:
                            shutil.copy(src_name, for_dir + sep + fle)
                        homeres[item]['add'].append((for_dir, fle, lib))
            # ~~ Parallel log files
            tasks = []
            mes = Messages(size=10, ncsize=ncsize)
            # ~~ Creates modules and objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if not homeres[item]['add']:
                print('        +> There is no need to compile any object')
            else:
                # ~~ Refactor some of the names
                if mod in cfg['ADDONES']:
                    print('        +> Refactoring the tree top ' + item)
                    refactor_sources(homeres[item]['add'],
                                     cmdf_files[mod][item], False)
                print('        +> Compile / Assemble / Link')
                ibar = 0
                pbar = \
                    ProgressBar(maxval=len(homeres[item]['add'])).start()
                for pth, obj, lib in homeres[item]['add']:
                    # path is where the .o and .mod will eventually
                    # be created, for each file
                    out = create_obj_files(\
                            obj, item,
                            {'libname': lib,
                             'type': get_principal_wrap_names(\
                                        pth + sep + obj)[0],
                             'path': pth},
                            mes, tasks,
                            bypass, homeres, verbose)
                    for _, _, err, _, file_name in out:
                        if err == '':
                            pbar.write('    - completed: ' + file_name, ibar)
                        else:
                            raise TelemacException(\
                              '\n        +> failed: ' + file_name + '\n' + err)
                        ibar = ibar + 1
                        pbar.update(ibar)
                # ~~> waiting for the remaining queued jobs to complete
                out = mes.flush_cmd(tasks)
                for _, _, err, _, file_name in out:
                    if err == '':
                        pbar.write('    - completed: ' + file_name, ibar)
                    else:
                        raise TelemacException(\
                              '\n        +> failed: ' + file_name + '\n' + err)
                    ibar = ibar + 1
                    pbar.update(ibar)
                pbar.finish()
                sleep(1)
            if item == 'HOMERE_MASCARET':
                compile_mascaret_dependencies(rebuild, homeres)
            # ~~ Creates libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            found_lib = True
            for lib in homeres[item]['deps']:
                prog = item.lower()
                if 'homere_' in item.lower():
                    prog = prog.split('homere_')[1]
                f = create_lib_files(lib, mod, item, prog,
                                     mes, tasks, bypass, homeres, verbose)
                # ~~> waiting for the remaining queued jobs to complete
                out = mes.flush_cmd(tasks)
                for _, _, err, _, file_name in out:
                    if err == '':
                        print('    - completed: ' + file_name)
                    else:
                        raise TelemacException(\
                             '\n        +> failed: ' + file_name + '\n')
                found_lib = found_lib and f
            if found_lib:
                print('        +> There is no need to package any library')
            else:
                sleep(1)
            # ~~ Creates executable ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if cfg['MODULES'][mod]['xexe'] != '':
                if mod not in cfg['ADDONES'].keys() and \
                        mod not in cfg['ODDONES'].keys():
                    found_exe = True
                    prog = item.lower()
                    if 'homere_' in prog:
                        prog = prog.split('homere_')[1]
                    found_exe = create_exe_files(item, mod,
                                                 prog, mes,
                                                 bypass, homeres, verbose)
                    if found_exe:
                        print('        +> There is no need to create '
                              'the associate executable')
            # ~~ Creates a python module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if cfg['MODULES'][mod]['xpyd'] != '':
                found_pyd = True
                f = ''
                for obj, lib in homeres[item]['add']:
                    if lib == mod + '.' + item:
                        f = obj
                for obj, lib in homeres[item]['tag']:
                    if lib == mod + '.' + item:
                        f = obj
                if f == '':
                    raise TelemacException(\
                          '\n        +> failed: ' + item
                          + ' not found in adds or tags of'
                          '\n                '
                          + repr(homeres[item]))
                prog = item.lower()
                if 'homere_' in prog:
                    prog = prog.split('homere_')[1]
                found_pyd = create_pyd_files(item, f, mod,
                                             prog, mes, bypass,
                                             homeres, verbose)
                if found_pyd:
                    print('        +> There is no need to create '
                          'the associate python module')
    # ~~ End of scans for all cmdf files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if not found:
        raise TelemacException(\
             'Could not find any cmdf file for config '
             + cfgname
             + '. You may have to use the --rescan option')
