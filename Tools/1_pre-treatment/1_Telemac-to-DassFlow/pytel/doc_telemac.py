#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium

   @brief Scripts to compile the telemac-mascaret documentation
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
from os import chdir, remove, walk, path, linesep, sep, mkdir, stat
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from string import Template
# ~~> dependencies towards other pytel/modules
from utils.messages import Messages
from utils.files import get_file_content
from utils.exceptions import TelemacException
from config import add_config_argument, update_config, CFGS

MODULE_LIST = ['artemis',
               'stbtel',
               'sisyphe',
               'postel3d',
               'telemac2d',
               'telemac3d',
               'tomawac',
               'waqtel',
               'telapy',
               'mascaret',
               'gaia',
               'nestor',
               'khione',
               'coupling']

MISC_LIST = ['developer_guide',
             'software_quality_plan',
             'TelemacDocTemplate',
             'git_guide',
             'doxydocs',
             'doxypydocs',
             'notebooks']

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#
__author__ = "Yoann Audouin"
__date__ = "$21-Sep-2012 16:51:09$"

def clean_doc(doc_dir, fullclean):
    """
    @brief Remove latex temporary files

    @param doc_dir Directory containing the main tex file
    @param fullclean If Yes will remove the pdf file as well
    """
    _, _, files = next(walk(doc_dir))
    for fle in files:
        if fle.endswith((".bbl", ".blg", ".aux", ".out", ".toc", ".log",
                         ".nlo", "~", "idx", "ptc")):
            remove(fle)
        if fullclean and fle.endswith(".pdf"):
            remove(fle)

def generate_list_variables(doc_dir):
    """
    Generates in latex format list of variables for each telapy module

    @param doc_dir (str) Path of the documentation
    """
    try:
        import _api
        api = sys.modules['_api']
        api_avail = True
    except ImportError:
        print("Api nor available generating empty list of variable")
        api_avail = False

    list_variable_template = Template(r"""
\begin{longtable}{|p{.40\textwidth} | p{.60\textwidth}|}
\hline
Variable name & Definition \tabularnewline
\hline
\hline
$content
\hline
\caption{Accessible variables through the API for $mod}
\end{longtable}""")

    for module in ['t2d', 't3d', 'sis', 'art', 'wac']:

        content = ''
        if api_avail:
            print(" ~> Getting var for", module)
            # Using direct acces to api (no need to instanciate)
            mod_handle_var = getattr(api, "api_handle_var_"+module)
            nb_var = getattr(mod_handle_var, "nb_var_"+module)
            var_len = getattr(mod_handle_var, module+"_var_len")
            info_len = getattr(mod_handle_var, module+"_info_len")
            get_var_info = getattr(mod_handle_var,
                                   "get_var_info_{}_d".format(module))
            set_var_list = getattr(mod_handle_var,
                                   "set_var_list_{}_d".format(module))

            # Building array of variable info
            error = set_var_list()
            if error != 0:
                raise TelemacException(
                    "Error while setting var list: {}".format(error))

            # Getting info for each variable
            for i in range(nb_var):
                tmp_varname, tmp_varinfo, error = \
                        get_var_info(i+1, var_len, info_len)
                if error != 0:
                    raise TelemacException(
                        "Error when getting info for var {}".format(i+1))
                varname = b''.join(tmp_varname).decode('utf-8').strip()
                varinfo = b''.join(tmp_varinfo).decode('utf-8').strip()

                # Adding new line the array
                content += '\n{} & {}\\tabularnewline'\
                           .format(varname.replace('_', r'\_'),
                                   varinfo.replace('_', r'\_'))

        file_name = path.join(doc_dir, 'latex',
                              'var_list_{}.tex'.format(module))
        with open(file_name, 'w') as ffile:
            ffile.write(list_variable_template.substitute(
                mod=module, content=content))

def compiletex(texfile, version, verbose):
    """
    @brief Full procedure for compiling a LaTeX file
             .i.e pdflatex,bibtex,pdflatex,pdflatex
    @param texfile Name of the main LaTex file
    @param version Version of the code/documentation
    @param verbose If yes display pdflatex listing
    """

    if verbose:
        tex_cmd = "pdflatex --jobname={tex}_{version} {tex}.tex"\
                  .format(tex=texfile, version=version)
        bib_cmd = "bibtex {}_{}.aux".format(texfile, version)
    else:
        tex_cmd = \
           "pdflatex --jobname={tex}_{version} {tex}.tex >latex_run.log 2>&1"\
                  .format(tex=texfile, version=version)
        bib_cmd = "bibtex {}_{}.aux >bib_run.log 2>&1".format(texfile, version)


    # First compilation
    mes = Messages(size=10)
    tail, code = mes.run_cmd(tex_cmd, False)

    if code != 0:
        if verbose:
            log = ''
        else:
            log = '\noutput:\n{}'.format(\
                    '\n'.join(get_file_content('latex_run.log')[-20:]))
        raise TelemacException(\
            'Latex compilation failed\n{}\n{}'\
            .format(tail, log))

    # Bibtex compilation
    tail, code = mes.run_cmd(bib_cmd, False)

    # Forcing bibtex to pass
    code = 0

    if code != 0:
        if verbose:
            log = ''
        else:
            log = '\noutput:\n{}'.format(\
                    '\n'.join(get_file_content('latex_run.log')[-20:]))
        raise TelemacException(\
            'Latex compilation failed\n{}\n{}'\
            .format(tail, log))

    # Second compilation
    tail, code = mes.run_cmd(tex_cmd, False)

    if code != 0:
        if verbose:
            log = ''
        else:
            log = '\noutput:\n{}'.format(\
                    '\n'.join(get_file_content('bib_run.log')[-20:]))
        raise TelemacException(\
            'Latex compilation failed\n{}\n{}'\
            .format(tail, log))

    # Third compilation
    tail, code = mes.run_cmd(tex_cmd, False)

    if code != 0:
        if verbose:
            log = ''
        else:
            log = '\noutput:\n{}'.format(\
                    '\n'.join(get_file_content('latex_run.log')[-20:]))
        raise TelemacException(\
            'Latex compilation failed\n{}\n{}'\
            .format(tail, log))

#
def create_case_list_file(doc_dir, val_dir, cfg_val, cleanup):
    """
    @brief Creates the CASELIST.tex which includes
          all the test cases tex file

    @param doc_dir Path to directory containing the main LaTeX file
    @param val_dir Path to directory containing validation test cases
    @param cfg_val list of path for the examples
    @param cleanup If yes clean up the temporay files instead
                   of creating the CASELIST.Tex file

    @return the list of cases that where missing the .tex file
    """
    case_list_file = path.join(doc_dir, 'latex', 'CASELIST.tex')
    # Creating latex folder if not there
    # often the case with git (empty folder are not created)
    if not path.exists(path.join(doc_dir, 'latex')):
        mkdir(path.join(doc_dir, 'latex'))
    skipped_cases = []
    if cleanup:
        if path.exists(case_list_file):
            remove(case_list_file)
    else:
        # Remove the file if it is already there
        if path.exists(case_list_file):
            remove(case_list_file)
        with open(case_list_file, 'w') as fobj:
            # Loop on all test cases
            for case in sorted(cfg_val):
                if not path.exists(path.join(val_dir, case, 'doc',
                                             case+".tex")):
                    skipped_cases.append(case)
                else:
                    txt = linesep + r'\subincludefrom{' + \
                        val_dir.replace(sep, '/') +\
                        '/' + case + '/' + 'doc' +\
                        '/' + '}{' + case + '}' + \
                        linesep + r'\clearpage' + linesep
                    fobj.write(txt)
    return skipped_cases


def generate_ref_from_dict(exe_path, dictionary, latex_file, lng,
                           cleanup, verbose):
    """
    @brief Generate the Latex file for the
            reference manual from the dictionary

    @param exe_path (string) Path to homere_damocles executable
    @param dictionary (string) Path to the dictionary to read
    @param latex_file (string) Name of the outpu latex file that will
                         contain the reference manual
    @param lng (int) Language for the reference manual
                 1: French
                 2: English
    @param cleanup (boolean) If True removing genrated files
    @param verbose (boolean) If True display command
    """
    # Building input parameter file
    latex_dir = path.dirname(latex_file)
    if not path.exists(latex_dir):
        mkdir(latex_dir)
    param_file = path.join(latex_dir, 'gen_ref.par')
    log_file = path.join(latex_dir, 'gen_ref.log')
    # Cleanup
    if cleanup:
        if path.exists(param_file):
            remove(param_file)
        if path.exists(log_file):
            remove(log_file)
        if path.exists(latex_file):
            remove(latex_file)
    else:
        # Creating parameter file for damocles
        with open(param_file, 'w') as f:
            f.write('LATEX'+'\n')
            f.write(dictionary+'\n')
            f.write(latex_file+'\n')
            f.write(lng+'\n')
        # Removing LaTeX file if already there
        if path.exists(latex_file):
            remove(latex_file)
        # Run Fortran program
        mes = Messages(size=10)
        if verbose:
            cmd = "{} < {}".format(exe_path, param_file)
        else:
            cmd = "{} < {} >{} 2>&1".format(exe_path, param_file, log_file)
        if verbose:
            print(cmd)
        _, code = mes.run_cmd(cmd, False)
        if code != 0:
            if verbose:
                log = ''
            else:
                log = '\n\nHere is the log:\n'\
                      + '\n'.join(get_file_content(log_file))
            raise TelemacException(\
                    'Could not generated data from dictionary '
                    + '{}'.format(log))

def compile_doc(doc_dir, doc_name, version, cleanup, fullcleanup, verbose):
    """
    @brief Compile the telemac-mascaret documentation

    @param doc_dir Directory containing the main LaTeX file
    @param doc_name Name of the main LaTeX file
    @param version Version of the code/documentation
    @param cleanup If yes remove temporary files
    @param fullcleanup If yes does cleanup + remove pdf
    @param verbose If yes display pdflatex listing
    """
    chdir(doc_dir)
    if cleanup or fullcleanup:
        clean_doc(doc_dir, fullcleanup)
        print('   - Cleaned up folder '+doc_dir+'\n')
    else:
        # removing pdflatex temporary files
        clean_doc(doc_dir, False)
        # compiling the texfile
        compiletex(doc_name, version, verbose)

def generate_notebook_html(doc_dir, notebook_dir, verbose):
    """
    Generate an html layout of the notebooks using ipython nbconvert
    Than coying back the file into doc_dir

    @param doc_dir (string) Path to the folder that will contain the html
                            version of the documentation
    @param notebook_dir (string) Path to the notebooks
    @param verbose (bool) If True more verbose
    """
    # Creating doc folder if necessary
    if not path.exists(doc_dir):
        mkdir(doc_dir)

    # Running convertion in notebook folder
    # Gathering all html files
    for root, subdirs, files in walk(notebook_dir):
        # Creating subfolders in advance
        for subdir in subdirs:
            if ".ipynb_checkpoint" in root:
                continue
            out_dir = path.join(doc_dir + root.replace(notebook_dir, ''),
                                subdir)
            if not path.exists(out_dir):
                mkdir(out_dir)
        for ffile in files:
            if ffile.endswith("ipynb"):
                # Skipping notebook tmp folders
                if ".ipynb_checkpoint" in root:
                    continue
                notebook = path.join(root, ffile)
                out_dir = doc_dir + root.replace(notebook_dir, '')
                if verbose:
                    log_lvl = 'DEBUG'
                else:
                    log_lvl = 'ERROR'
                cmd = "jupyter nbconvert --to html --log-level={log_lvl} "\
                      "--output-dir={out_dir} --output={output} {nb}"\
                       .format(log_lvl=log_lvl, out_dir=out_dir,
                               output="tmp.html", nb=notebook)
                print("   ~> Converting "+\
                        path.join(root.replace(notebook_dir, '')[1:], ffile))
                if verbose:
                    print(cmd)
                # Running convertion
                mes = Messages(size=10)
                tail, code = mes.run_cmd(cmd, bypass=False)

                if code != 0:
                    raise TelemacException('nbconvert failed\n {}'.format(tail))

                tmp_file = path.join(out_dir, 'tmp.html')
                out_file = path.join(out_dir, ffile[:-5] + "html")

                # Replacing .ipynb in content of file by .html
                with open(tmp_file, 'r') as f:
                    content = f.read()

                remove(tmp_file)
                with open(out_file, 'w') as f:
                    f.write(content.replace(".ipynb", ".html"))


def compile_doxygen(doxy_file, verbose):
    """
    Compile a doxygen documentation

    @param doxy_file name of the doxygen file to use
    @param verbose If True display doxygen listing
    """

    doxy_dir = path.dirname(doxy_file)

    chdir(doxy_dir)
    if verbose:
        cmd = "doxygen {}".format(doxy_file)
    else:
        cmd = "doxygen {} >Doxygen_run.log 2>&1 ".format(doxy_file)

    print("   ~> Generating doxygen documentation for " + doxy_file)
    # Running convertion
    mes = Messages(size=10)
    tail, code = mes.run_cmd(cmd, bypass=False)

    if code != 0:
        if verbose:
            log = ''
        else:
            log = '\noutput:\n{}'.format(\
                    '\n'.join(get_file_content('Doxygen_run.log')[-20:]))
        raise TelemacException(\
            'Doxygen failed\n{}\n{}'\
            .format(tail, log))

    html_file = path.join(doxy_dir, 'html', 'index.html')

    print("   To see documentation run (replace firefox by your "\
          "internet browser):\n   firefox {}".format(html_file))

    doxy_warning_log = path.join(doxy_dir, 'Doxygen_warning.log')

    if stat(doxy_warning_log).st_size != 0:
        print('There seems to be some doxygen warnings see:\n{}'\
              .format(doxy_warning_log))

def generate_doxygen(doxydoc, verbose):
    """
    Generate the Doxygen documentation (In html) for the Python and the sources

    @param doxydoc name of the doxygen folder to use (doxydocs or doxypydocs)
    @param verbose If True display doxygen listing
    """

    # Compiling sources doxygen
    doxy_dir = path.join(CFGS.get_root(),
                         'documentation',
                         doxydoc)
    doxy_file = path.join(doxy_dir, 'Doxyfile')

    if doxydoc == "doxypydocs":
        try:
            import doxypypy
        except ImportError:
            raise TelemacException("doxypypy is mandatory to compile doxygen "
                                   "for Telemac scripts")

    compile_doxygen(doxy_file, verbose)

def main():
    """
    Main program for the compilation of the documentation of
    the telemac-mascaret system
    """
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nLoading Options and Configurations\n'+'~'*72+'\n')
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        description=('''\n
By Default all the documentation are generated\n
use the options --validation/reference/user/release/theory to compile only one
        '''))
    parser = add_config_argument(parser)
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        dest="verbose", default=False,
        help="Will display listing for all commands")
    parser.add_argument(
        "-m", "--modules",
        dest="modules", default='',
        help="specify the list modules (, separated), default is all of them "+
             "from {"+",".join(MODULE_LIST)+"}")
    parser.add_argument(
        "-M", "--misc",
        dest="misc", default='',
        help="specify the list of misc documentation (, separated) to compile, "
             "default is all of them "+
             "from {"+",".join(MISC_LIST)+"}")
    parser.add_argument(
        "--validation", action="store_true",
        dest="validation", default=False,
        help="Will generate the validation documentation")
    parser.add_argument(
        "--case-list",
        dest="case_list", default='',
        help="List of cas to include in the validation documentation"
             "separated by ',' (default all of them)")
    parser.add_argument(
        "--reference", action="store_true",
        dest="reference", default=False,
        help="Will generate the reference documentation")
    parser.add_argument(
        "--user", action="store_true",
        dest="user", default=False,
        help="Will generate the user documentation")
    parser.add_argument(
        "--release", action="store_true",
        dest="release_note", default=False,
        help="Will generate the release note")
    parser.add_argument(
        "--theory", action="store_true",
        dest="theory_guide", default=False,
        help="Will generate the theory guide")
    parser.add_argument(
        "--clean", action="store_true",
        dest="cleanup", default=False,
        help="Will remove all temporary file generated by pdflatex")
    parser.add_argument(
        "--fullclean", action="store_true",
        dest="fullcleanup", default=False,
        help="Same as clean but removes the pdf as well")

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Environment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    options = parser.parse_args()
    update_config(options)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Works for all configurations unless specified ~~~~~~~~~~~~~~~
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Compile the valiation documentation
    doall = not (options.validation or options.user or options.reference
                 or options.release_note or options.theory_guide)
    cfg = CFGS.configs[CFGS.cfgname]
    # still in lower case
    root = CFGS.get_root()
    # Get what i to be compiled
    # By default everything if something is defined compiling only that
    if options.modules != '':
        module_list = options.modules.split(',')
        for module in module_list:
            if module not in MODULE_LIST:
                raise TelemacException("{} is not in list of modules ({})".format(module, ",".join(MODULE_LIST)))
    else:
        # all modules
        module_list = MODULE_LIST
    if options.misc != '':
        misc_list = options.misc.split(',')
        for misc in misc_list:
            if misc not in MISC_LIST:
                raise TelemacException("{} is not in list of misc ({})".format(misc, ",".join(MISC_LIST)))
        module_list = []
    else:
        # all docs
        misc_list = MISC_LIST
        # If a module was specified or a specific documentation for modules
        # not compiling Misc documentation
        if options.modules != '' or not doall:
            misc_list = []

    CFGS.compute_vnv_info()

    # Get version in config if it exist use trunk otherwise
    version = cfg.get('version', 'trunk')

    # Initialise output message
    output_mess = '\n\n'
    # Look on all the modules for the documentation
    for code_name in module_list:
        print('\nCompilation of the documentation for ' + code_name
              + '\n'+'~'*72)
        # list of what to do for the module
        todo = []
        if options.validation or doall:
            if code_name not in ['artemis', 'telapy', 'coupling']:
                # Building Validation LaTeX file
                doc_dir = path.join(root, 'documentation',
                                    code_name, 'validation')
                chdir(doc_dir)
                if options.case_list != '':
                    list_of_case = options.case_list.split(',')
                else:
                    list_of_case = list(cfg['VALIDATION'][code_name].keys())
                    list_of_case.remove('path')
                skiped_case = \
                    create_case_list_file(
                        doc_dir,
                        cfg['VALIDATION'][code_name]['path'],
                        list_of_case,
                        options.cleanup or options.fullcleanup)
                for case in skiped_case:
                    output_mess += r'   - /!\ Missing LaTeX file for ' + \
                                   case+'\n'
                todo.append('validation')
        if options.reference or doall:
            if code_name not in ['telapy', 'mascaret', 'nestor', 'coupling']:
                # Path to the dictionary
                dictionary = path.join(root, 'sources', code_name,
                                       code_name+'.dico')
                # Path to latex File
                latex_file = path.join(root, 'documentation',
                                       code_name, 'reference',
                                       'latex', 'Corpus.tex')
                # English only for now
                lng = '2'
                # Path to bin directory
                exe_path = path.join(\
                        root, 'builds', CFGS.cfgname,
                        'bin', 'damocles'+cfg['SYSTEM']['sfx_exe'])
                generate_ref_from_dict(\
                        exe_path, dictionary, latex_file, lng,
                        options.cleanup or options.fullcleanup,
                        options.verbose)
                todo.append('reference')
        if options.user or doall:
            if code_name not in ['mascaret']:
                # Normal Compilation of a LaTeX file
                todo.append('user')
        if options.theory_guide or doall:
            # theory guide only available for telemac3d
            if code_name in ['telemac3d', 'mascaret', 'waqtel']:
                todo.append('theory_guide')
        for doc_type in todo:
            print('\n     ~> Compilation of the {} documentation'\
                  .format(doc_type))
            doc_dir = path.join(root, 'documentation',
                                code_name, doc_type)
            chdir(doc_dir)
            # Check if the file exist
            if path.exists(path.join(doc_dir,
                                     code_name + "_" + doc_type + ".tex")):
                if code_name == 'telapy' and \
                   not (options.cleanup or options.fullcleanup):
                    # Running small script to generate list of variables for api
                    generate_list_variables(doc_dir)
                compile_doc(doc_dir, code_name+'_'+doc_type,
                            version,
                            options.cleanup, options.fullcleanup,
                            options.verbose)
            else:
                raise TelemacException(\
                        "   - Error for {} {}, {}.tex "
                        "not found ".format(code_name,
                                            path.basename(doc_dir),
                                            code_name+"_"+doc_type))
            if not (options.cleanup or options.fullcleanup):
                output_mess += '   - Created %s_%s_%s.pdf\n' % \
                              (code_name, doc_type, version)
    # List of the other documentation
    print('\nCompilation of the documentation for Misc'
          + '\n'+'~'*72)
    for doc in misc_list:
        print('\n     ~> Compilation of the {} documentation'.format(doc))
        doc_dir = path.join(root, 'documentation',
                            'Misc', doc)

        if doc == 'notebooks':
            notebook_dir = path.join(root, 'notebooks')
            doc_dir = path.join(root, 'documentation', doc)
            if not (options.fullcleanup or options.cleanup):
                generate_notebook_html(doc_dir, notebook_dir, options.verbose)
        elif doc in ['doxydocs', 'doxypydocs']:
            if not (options.fullcleanup or options.cleanup):
                generate_doxygen(doc, options.verbose)
        else:
            chdir(doc_dir)
            if path.exists(path.join(doc_dir, doc + ".tex")):
                compile_doc(doc_dir, doc,
                            version,
                            options.cleanup, options.fullcleanup,
                            options.verbose)
            else:
                raise TelemacException(\
                        "   - Error in {}, {}.tex "
                        "not found ".format(path.basename(doc_dir), doc))

        if not (options.cleanup or options.fullcleanup):
            if doc not in ['notebooks', 'doxydocs', 'doxypydocs']:
                output_mess += '   - Created %s for %s.pdf\n' % \
                              (doc, version)
            else:
                output_mess += '   - Created %s_%s.pdf\n' % \
                              (doc, version)


    print(output_mess)
    print('\n\n'+'~'*72)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)


if __name__ == "__main__":
    main()
