"""@author TELEMAC-MASCARET Consortium

   @brief Run a converions of mesh files using stbtel
"""
from os import path
import re
from utils.exceptions import TelemacException

CHECK_RES_TEMPLATE =\
"""
        # {title}
        self.check_epsilons('{file1}',
                            '{file2}',
                            eps={eps})
"""

GET_FILE_TEMPLATE =\
"""        {var} = self.get_study_file('{target}')
        res_{var} = TelemacFile({var})
"""

PLOT1D_HISTORY_TEMPLATE =\
"""
        #Plotting {variable} on {points} over records {records}
        vnv_plot1d_history(\\
                '{variable}',
                res_{file_var},
                '{variable}',
                points=[{points}],
                fig_size={fig_size},
                fig_name='{fig_name}')

"""

PLOT1D_TIMESERIE_POLY_TEMPLATE =\
"""
        # Plotting {variable} over polyline over records {records}
        vnv_plot1d_polylines(\\
                '{variable}',
                res_{file_var},
                poly=[{poly}],
                record={records},
                fig_size={fig_size},
                fig_name='{fig_name}')
"""

PLOT_VERTICAL_TEMPLATE =\
"""
        # Plotting vertical split
        vnv_plot2d('{variable}',
                   res_{file_var},
                   poly=[{poly}],
                   record={record},
                   filled_contours=True,
                   fig_size={fig_size},
                   fig_name='{fig_name}')
"""

PLOT_HORIZONTAL_TEMPLATE =\
"""
        # Plotting horizontal split
        vnv_plot2d('{variable}',
                   res_{file_var},
                   plane={plane},
                   record={record},
                   filled_contours=True,
                   fig_size={fig_size},
                   fig_name='{fig_name}')
"""


PLOT2D_MESH_TEMPLATE =\
"""
        #Plotting mesh
        vnv_plot2d('{variable}',
                   res_{file_var},
                   plot_mesh=True,
                   fig_size={fig_size},
                   fig_name='{fig_name}')

"""
PLOT2D_VAR_TEMPLATE =\
"""
        # Plotting {variable} at {record}
        vnv_plot2d('{variable}',
                   res_{file_var},
                   record={record},
                   fig_size={fig_size},
                   fig_name='{fig_name}')

"""

PLOT2D_3D_SCALAR_TEMPLATE =\
"""
        # Plotting 3d scalar map for {variable} at {record}
        plot_3d_scalar_map({file_var},
                           '{variable}',
                           record={record},
                           fig_size={fig_size},
                           fig_name='{fig_name}')

"""

ADD_STUDY_TEMPLATE =\
"""
        # {title}
        self.add_study('{name}',
                       '{module}',
                       '{cas}')

"""

ADD_CMD_TEMPLATE =\
"""
        # {title}
        self.add_command('{name}',
                         '{command}')

"""

ADD_PARA_STUDY_TEMPLATE =\
"""
        # {title}
        cas = TelemacCas('{cas_seq}', get_dico('{module}'))
        cas.set('PARALLEL PROCESSORS', {ncsize}){add_set}

        self.add_study('{name}',
                       '{module}',
                       '{cas}',
                       cas=cas)

        del cas

"""

PY_TEMPLATE =\
'''
"""
Validation script for gouttedo
"""
from vvytel.vnv_study import AbstractVnvStudy
from execution.telemac_cas import TelemacCas, get_dico

class VnvStudy(AbstractVnvStudy):
    """
    Class for validation of Thompson
    """

    def _init(self):
        """
        Defines the general parameter
        """
        self.rank = {rank}
        self.tags = [{tags}]

    def _pre(self):
        """
        Defining the studies
        """
        {pre}

    def _check_results(self):
        """
        Post-treatment processes
        """
        {check_results}

    def _post(self):
        """
        Post-treatment processes
        """
        {post}'''

def xml2py_parser(subparser):
    """
    Add arguments for a xml2py conversion

    @param subparser (argumentParser) argument parser

    @return (argumentParser) the updated argument parser
    """
    parser = subparser.add_parser('xml2py',\
            help='Convert a xml for validation into a Python for validation')
    parser.add_argument(
        "input_file", default="",
        help="name of the input file also defines the input format")
    parser.add_argument(
        "--skip", dest="skip", action='store_true', default=False,
        help="If True will skip unhandled features but still do the conversion")
    parser.add_argument(
        "-o",
        dest="out_file", default=None,
        help="Force the name of the output file default will be vnv_{xref}.py"\
             "where xref is the attribute of xml")

    return subparser

def records2pylist(string, file_var):
    """
    Will convert a string that looks like:
      - [i] were i can be -1 and must be replace by file_var.ntimestep
      - [i:j] will become range(i, j)
      - [i:j:] will become range(i, j, k)

    @param string (str) String ro parse
    @param file_var (str) Name of the python variable to use

    @returns (str) The new string
    """
    # Special case (only in gouttedo.xml...)
    if string == "[0:0]":
        return '[0]'
    string2 = string[1:-1].replace('-1', 'res_{}.ntimestep'.format(file_var))
    if ':' in string2:
        if string2.count(':') == 1:
            i, j = string2.split(':')
            return 'range({}, {})'.format(i, j)
        else:
            i, j, k = string2.split(':')
            return 'range({}, {}, {})'.format(i, j, k)
    elif ';' in string2:
        records = string2.split(';')
        return '[{}]'.format(', '.join(records))
    else:
        return string2

def extract2pylist(string):
    """
    Convert an extraction string (from xml) into its Python equivalent

    @param string (str) the string to convert

    @returns (str)  The converted string
    """
    poly = string.replace(')(', '), (')\
                 .replace(';', ', ')\
                 .replace(':', '')\
                 .replace('(', '[')\
                 .replace(')', ']')
    return poly

def xml2py(xml_file, out_file=None, skip=False):
    """
    Converter a validation xml file into a Python validation file

    @param xml_file (str) Path to the xml file to convert
    @param out_file (str) If given it will be name of the converter file
    otherwise it will be named vnv_{name}.py where name is the xref given with
    the xml
    @param skip (bool) If True will not crashed when unhandled tags are found but will just skip them instead
    """

    import xml.etree.ElementTree as XML
    with open(xml_file, 'r') as f:
        xml_tree = XML.parse(f)  # may need to try first and report error
        xml_root = xml_tree.getroot()

    rank = xml_root.attrib['rank']
    tags = xml_root.attrib['tags']
    name = xml_root.attrib['xref']
    tags_str = ", ".join(["'{}'".format(tag) for tag in tags.split(';')])
    pre_str = ''
    post_str = ''
    check_results_str = ''

    if out_file is None:
        name = path.splitext(path.basename(xml_file))[0]
        py_file = path.join(path.dirname(xml_file), 'vnv_'+name+'.py')
    else:
        py_file = out_file

    # USed to store cast variables
    var = {}
    # Used to store files used for post-treatment this will containe the name
    # of the Python variable associated with the file
    files = {}
    # Looping on all xml elements
    for xml_child in xml_root:
        # Action xml element
        if 'action' in xml_child.tag:
            # code action correspond to doing a add_study in Python
            if 'code' in xml_child.attrib:
                if xml_child.attrib['code'] == 'exec':
                    xref = xml_child.attrib['xref']
                    do = xml_child.attrib['do']
                    title = xml_child.attrib['title']
                    pre_str += ADD_CMD_TEMPLATE.format(\
                            title=title,
                            name=xref,
                            command=do)
                else:
                    # If ncsize is presetn modifying steering to add parallel
                    # processors keyword
                    if 'ncsize' in xml_child.attrib:
                        # Parallel file to create
                        cas_file = xml_child.attrib['target']
                        module = xml_child.attrib['code']
                        title = xml_child.attrib['title']
                        ncsize = xml_child.attrib['ncsize']
                        xref = xml_child.attrib['xref']

                        root, ext = path.splitext(cas_file)
                        cas_para = root + "_par"  + ext
                        if module == 'artemis' and 'set' in xml_child.attrib:
                            # Setting solver to mumps for artemis
                            add_set = "\n        cas.set('SOLVER', 9)\n"
                        else:
                            add_set = ""
                        pre_str += ADD_PARA_STUDY_TEMPLATE.format(\
                                name=xref,
                                title=title,
                                module=module,
                                cas=cas_para,
                                cas_seq=cas_file,
                                ncsize=ncsize,
                                add_set=add_set)
                    else:
                        # Adding study
                        cas_file = xml_child.attrib['target']
                        module = xml_child.attrib['code']
                        title = xml_child.attrib['title']
                        xref = xml_child.attrib['xref']
                        pre_str += ADD_STUDY_TEMPLATE.format(\
                                name=xref,
                                title=title,
                                module=module,
                                cas=cas_file)
        # Cast are mainly used for checking epsilons
        # We need to extract the variables cast (usually from a results file)

        elif 'cast' in xml_child.tag:
            for xml_cast in xml_child:
                # Return element containe the call the checkval and the epsilons
                if xml_cast.tag == 'return':
                    # Getting the title and the epsilons used for the diff
                    title = xml_cast.attrib['title']
                    proc = re.match(r"checkval\([\w\-]+,\s*(?P<eps>\[[0-9.,eE\- ]*\])(,\s*norm='[\w]*')?\)", xml_cast.attrib['fail'])
                    eps = proc.group('eps').replace(' ', '').replace(',', ', ')
                # Python element are not handled
                elif xml_cast.tag == 'python':
                    if not skip:
                        raise TelemacException(\
                                "Python tag not handled by converter")
                    else:
                        print("  ~> Skipping (Python tag): ",
                              xml_child.attrib['xref'])
                        continue
                # If the element contains a target element this is casting data
                # do we extract what we need
                elif 'target' in xml_cast.attrib:
                    # Add to vars (this will map variable defines into the xml
                    # name to their value)
                    name = xml_cast.tag
                    target = xml_cast.attrib['target']
                    var[name] = target
                else:
                    # Should be the mapdiff one that contains the name of the
                    # two variable we are comparing
                    proc = re.match(r'mapdiff\((?P<file1>[\w\-]+),\s*(?P<file2>[\w\-+]+)(,\s*relative=False|,\s*relative=True|,\s*noname=True)*\)', xml_cast.attrib['vars'])
                    file1 = proc.group('file1')
                    file2 = proc.group('file2')
            # The loop before should have got us all we need
            check_results_str += CHECK_RES_TEMPLATE.format(\
                                     title=title,
                                     file1=var[file1],
                                     file2=var[file2],
                                     eps=eps)

        # 2D plot xml element
        elif 'plot2d' in xml_child.tag:
            fig_name = xml_child.attrib['xref']
            # fig size is often given but not always
            if 'size' in xml_child.attrib:
                fig_size = xml_child.attrib['size'].replace(';', ', ')
            else:
                fig_size = None
            split_type = xml_child.attrib.get('type', '')
            # Extract and time can be given in either the main element or in
            # the layers
            record = xml_child.attrib.get('time', '[0]')
            extract = None
            if 'extract' in xml_child.attrib:
                extract = xml_child.attrib['extract']

            if len(xml_child) > 1:
                if not skip:
                    raise TelemacException("Not handling multiple layer")
                else:
                    print("  ~> Skipping (Multiple layers): ",
                          xml_child.attrib['xref'])
                    continue
            for xml_layer in xml_child:
                variable, plot_type = xml_layer.attrib['vars'].split(':')
                target = xml_layer.attrib['target']
                # Overwritting time as in can be in xml_child and layer
                if 'time' in xml_layer.attrib:
                    record = xml_layer.attrib.get('time')
                if 'extract' in xml_layer.attrib:
                    extract = xml_layer.attrib['extract']

                # If we have a new extraction source
                if target not in files:
                    # xref:TEGEO type
                    if ':' in target:
                        files[target] = target.replace(':', '_')\
                                              .replace('-', '_')\
                                              .lower()
                    # Directly a file name
                    else:
                        files[target] = "'{}'".format(target)

            if plot_type == 'map3d':
                if split_type == 'p-section':
                    post_str += PLOT2D_3D_SCALAR_TEMPLATE.format(\
                            file_var=files[target],
                            fig_name=fig_name,
                            fig_size=fig_size,
                            record=record[1:-1],
                            variable=variable.upper().strip())
            elif plot_type == 'mesh':
                if split_type == '':
                    post_str += PLOT2D_MESH_TEMPLATE.format\
                            (file_var=files[target],
                             variable=variable.upper().strip(),
                             fig_name=fig_name,
                             fig_size=fig_size)
            elif plot_type == 'map':
                if split_type == '':
                    post_str += PLOT2D_VAR_TEMPLATE.format(\
                             file_var=files[target],
                             fig_name=fig_name,
                             record=record[1:-1],
                             variable=variable.upper().strip(),
                             fig_size=fig_size)
                elif split_type == 'v-section':
                    post_str += PLOT_VERTICAL_TEMPLATE.format(\
                              file_var=files[target],
                              variable=variable.upper().strip(),
                              poly=extract2pylist(extract),
                              record=record[1:-1],
                              fig_name=fig_name,
                              fig_size=fig_size)
                elif split_type == 'p-section':
                    if extract is None:
                        plane = "0"
                    else:
                        # Plane as the format []#[i]
                        plane = extract.replace('[', '')\
                                       .replace(']', '')\
                                       .replace('#', '')

                    post_str += PLOT_HORIZONTAL_TEMPLATE.format(\
                              file_var=files[target],
                              variable=variable.upper().strip(),
                              plane=plane,
                              record=record[1:-1],
                              fig_name=fig_name,
                              fig_size=fig_size)
                else:
                    raise TelemacException("Unknow type: ", split_type)

        # Making a 1d plot
        elif 'plot1d' in xml_child.tag:
            # Converting a 1d plot
            fig_name = xml_child.attrib['xref']
            if 'size' in xml_child.attrib:
                fig_size = xml_child.attrib['size'].replace(';', ', ')
            else:
                fig_size = None
            split_type = xml_child.attrib.get('type', '')
            # extract and time can be either in the plot1d xml element or in
            # the layer xml element
            records = xml_child.attrib.get('time', '[0]')
            if 'extract' in xml_child.attrib:
                extract = xml_child.attrib['extract']

            if len(xml_child) > 1:
                if not skip:
                    raise TelemacException("Not handling multiple layers")
                else:
                    print("  ~> Skipping (Mulitple layers):",
                          xml_child.attrib['xref'])
                    continue
            for xml_layer in xml_child:
                variable, plot_type = xml_layer.attrib['vars'].split(':')
                target = xml_layer.attrib['target']
                # Overwritting time and extract
                if 'time' in xml_layer.attrib:
                    records = xml_layer.attrib['time']
                if 'extract' in xml_layer.attrib:
                    extract = xml_layer.attrib['extract']

                if target not in files:
                    # Building name of the python variable for the file name
                    files[target] = target.replace(':', '_')\
                                          .replace('-', '_')\
                                          .lower()

            if plot_type == 'line':
                if split_type == 'v-section':
                    records_str = records2pylist(records, files[target])

                    post_str += PLOT1D_TIMESERIE_POLY_TEMPLATE.format(\
                                     file_var=files[target],
                                     fig_name=fig_name,
                                     fig_size=fig_size,
                                     variable=variable.upper().strip(),
                                     records=records_str,
                                     poly=extract2pylist(extract))
                elif split_type == 'history':
                    records_str = records2pylist(records, files[target])
                    # TODO: Not handling 1d plot on a specific plane in 3d
                    if '#' in extract:
                        if not skip:
                            raise TelemacException(\
                                    'Unhandled plot1 history on 3d data')
                        else:
                            print("~> Skipping (plot1d history on 3d data): ",
                                  xml_child.attrib['xref'])
                            continue
                    post_str += PLOT1D_HISTORY_TEMPLATE.format(\
                                     file_var=files[target],
                                     fig_name=fig_name,
                                     fig_size=fig_size,
                                     variable=variable.upper().strip(),
                                     records=records_str,
                                     points=extract2pylist(extract))
                else:
                    if not skip:
                        raise TelemacException('Unhandled type for plot1d: ', split_type)
                    else:
                        print("~> Skipping (plot type): ",
                              xml_child.attrib['xref'])
                        continue

        elif xml_child.tag in ['deco', 'meta']:
            # Ignoring deco and meta for now
            pass
        else:
            raise TelemacException("Unhandled xml tag: ", xml_child.tag)

    # Adding opening of files in post
    if files != {}:
        files_str = "        # Getting files\n"
        for target, var in files.items():
            if ':' in target:
                files_str += GET_FILE_TEMPLATE.format(\
                        var=var,
                        target=target)
        post_str = files_str + post_str

    # Writing the file
    with open(py_file, 'w') as f:
        f.write(PY_TEMPLATE.format(\
                pre=pre_str, \
                post=post_str, \
                check_results=check_results_str,
                rank=rank,
                tags=tags_str))
