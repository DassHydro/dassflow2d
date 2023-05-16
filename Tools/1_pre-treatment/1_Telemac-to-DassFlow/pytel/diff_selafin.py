#!/usr/bin/env python3
"""
    @author TELEMAC-MASCARET Consortium
    @brief Reporting on differences between two SELAFIN files

"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
from os import path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import numpy as np
# ~~> dependencies towards other modules
from pretel.scan_selafin import ScanSelafin
from utils.exceptions import TelemacException
from compilation.parser_fortran import clean_quotes

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#

# _____                    _________________________________________
# ____/ Secondary Classes /________________________________________/
#

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

def main():
    """
    Main function of diffSELAFIN
    """
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nLoading Options and Configurations\n'+72*'~'+'\n')
    parser = ArgumentParser(\
        formatter_class=RawDescriptionHelpFormatter,
        description=('''\n
Reporting on differences between two SELAFIN files
        '''),
        usage=' (--help for help)\n---------\n      =>  '\
                '%(prog)s [options] file1.slf file2.slf\n---------')
    # ~~> Uselessly set to True as default ... may change in the future
    # ~~> The real switches
    parser.add_argument(\
        "--head", action="store_true",
        dest="head", default=False,
        help="Will print a statiscal differences between two SELARING files")
    parser.add_argument(\
        "--core", action="store_true",
        dest="core", default=False,
        help="Will print a statiscal differences between two SELARING files")
    parser.add_argument(\
        "--scan", action="store_true",
        dest="scan", default=False,
        help="Will print an individual summary for each file")
    parser.add_argument(\
        "-v", "--vars",
        dest="xvars", default=None,
        help= \
        "specify which variables should be differentiated (':'-delimited)")
    parser.add_argument(\
        "-f", "--from",
        dest="tfrom", default="1",
        help="specify the first frame included in the differentiation")
    parser.add_argument(\
        "-s", "--stop",
        dest="tstop", default="-1",
        help="specify the last frame included (negative from the end) "\
              "in the differentiation")
    parser.add_argument(\
        "-diff", "--step",
        dest="tstep", default="1",
        help="specify the step for the extraction of frames for "\
              "the differentiation")
    parser.add_argument(\
        "-e", "--epsilon",
        dest="epsilon", default="0",
        help="specify the threshold for which values are assumed the same")
    parser.add_argument(\
        "-b", "--bypass", action="store_true",
        dest="bypass", default=False,
        help="Will bypass certain mismatches between files")
    parser.add_argument(\
        "args", metavar='file1,file2',
        default='', nargs=2,
        help="operation: ( files1 - file2 )")
    options = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Double checks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    slf_file1 = options.args[0]
    if not path.exists(slf_file1):
        raise TelemacException('\nCould not find the file named: {}'.format(slf_file1))
    slf_file2 = options.args[1]
    if not path.exists(slf_file2):
        raise TelemacException('\nCould not find the file named: {}'.format(slf_file2))

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Initial scan ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    vrbls = options.xvars
    if options.xvars != None:
        vrbls = clean_quotes(options.xvars.replace('_', ' '))
    times = (int(options.tfrom), int(options.tstep), int(options.tstop))
    slf1 = ScanSelafin(slf_file1, times=times, vrs=vrbls)
    slf2 = ScanSelafin(slf_file2, times=times, vrs=vrbls)

    if options.scan:
        print('\n\nFirst file: '+slf_file1+'\n'+72*'~'+'\n')
        slf1.print_header()
        slf1.print_time_summary()
        print('\n\nSecond file: '+slf_file2+'\n'+72*'~'+'\n')
        slf2.print_header()
        slf2.print_time_summary()

    comparable = True

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Header differences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if options.head:
        print('\n\nHeader differences: \n'+72*'~'+'\n')

        # ~~> File formats
        if slf1.slf.file['endian'] != slf2.slf.file['endian']:
            print('\n  <> File ENDIANs:\n')
            if slf1.slf.file['endian'] == ">":
                print('     + '+slf_file1+' is BIG ENDIAN')
            else:
                print('     + '+slf_file1+' is LITTLE ENDIAN')
            if slf2.slf.file['endian'] == ">":
                print('     + '+slf_file2+' is BIG ENDIAN')
            else:
                print('     + '+slf_file2+' is LITTLE ENDIAN')
        if slf1.slf.file['float'] != slf2.slf.file['float']:
            print('\n  <> File FLOATs:\n')
            if slf1.slf.file['float'] == ('diff', 8):
                print('     + '+slf_file1+' is DOUBLE PRECISION')
            else:
                print('     + '+slf_file1+' is SINGLE PRECISION')
            if slf2.slf.file['float'] == ('diff', 8):
                print('     + '+slf_file2+' is DOUBLE PRECISION')
            else:
                print('     + '+slf_file2+' is SINGLE PRECISION')

    # ~~> File contents
    mes = '\n  <> List of variable names:\n'
    found = False
    cmn_vars = []
    mes = mes + '\n     + '+slf_file1
    for ivar in range(len(slf1.slf.varnames)):
        if slf1.slf.varnames[ivar] in slf2.slf.varnames:
            mes = mes + '\n        = '+slf1.slf.varnames[ivar]
            cmn_vars.append(slf1.slf.varnames[ivar])
        else:
            mes = mes + '\n        * '+slf1.slf.varnames[ivar]
            found = True
    mes = mes + '\n     + '+slf_file2
    for ivar in range(len(slf2.slf.varnames)):
        if slf2.slf.varnames[ivar] in slf1.slf.varnames:
            mes = mes + '\n        = '+slf2.slf.varnames[ivar]
        else:
            mes = mes + '\n        * '+slf2.slf.varnames[ivar]
            found = True
    if found and options.head:
        print(mes)
    if not cmn_vars:
        comparable = False
        print('\n  /!\\ no common variables. The files are not comparables.\n')

    # ~~> File reference dates and times
    if options.head:
        if max(np.array(slf1.slf.datetime) - np.array(slf2.slf.datetime)) > 0:
            print('\n  <> Different reference dates:')
            print('     + '+slf_file1+': '+repr(slf1.slf.datetime))
            print('     + '+slf_file2+': '+repr(slf2.slf.datetime))

    # ~~> File time frames
    mes = '\n  <> List of time frames:\n'
    found = False
    times0 = []
    times1 = []
    times2 = []
    # ~~> check if sorted times
    it1 = 1
    if len(slf1.slf.tags['times']) > 1:
        for it1 in range(len(slf1.slf.tags['times']))[1:]:
            if slf1.slf.tags['times'][it1] <= slf1.slf.tags['times'][it1-1]:
                break
        if slf1.slf.tags['times'][it1] > slf1.slf.tags['times'][it1-1]:
            it1 += 1
    it2 = 1
    if len(slf2.slf.tags['times']) > 1:
        for it2 in range(len(slf2.slf.tags['times']))[1:]:
            if slf2.slf.tags['times'][it2] <= slf2.slf.tags['times'][it2-1]:
                break
        if slf2.slf.tags['times'][it2] > slf2.slf.tags['times'][it2-1]:
            it2 += 1
    # ~~> correct if not bypassed
    if options.bypass and \
        len(slf1.slf.tags['times']) == len(slf2.slf.tags['times']):
        times0 = list(range(len(slf1.slf.tags['times'])))
    else:
        diff = np.setdiff1d(slf1.slf.tags['times'][:it1],
                            slf2.slf.tags['times'][:it2])
        if diff.size != 0:
            found = True
            mes = mes + '\n     + frames only in '+slf_file1+' : '+\
                    ', '.join(['{0:.2f}'.format(i) for i in diff])
        diff = np.setdiff1d(slf2.slf.tags['times'][:it2],
                            slf1.slf.tags['times'][:it1])
        if diff.size != 0:
            found = True
            mes = mes + '\n     + frames only in '+slf_file2+' : '+\
                    ', '.join(['{0:.2f}'.format(i) for i in diff])
        diff = np.intersect1d(slf1.slf.tags['times'][:it1],
                              slf2.slf.tags['times'][:it2])
        if diff.size != 0:
            mes = mes + '\n     + frames in both files: '+\
                  ', '.join([str(i) for i in diff])
            times1 = np.searchsorted(slf1.slf.tags['times'][:it1], diff)
            slf1.slf.tags['times'] = slf1.slf.tags['times'][times1]
            for time in range(len(slf1.slf.tags['cores']))[it1:]:
                slf1.slf.tags['cores'].remove(slf1.slf.tags['cores'][-1])
            for time in range(len(slf1.slf.tags['cores']))[::-1]:
                if time not in times1:
                    slf1.slf.tags['cores'].remove(slf1.slf.tags['cores'][time])
            times2 = np.searchsorted(slf2.slf.tags['times'][:it2], diff)
            slf2.slf.tags['times'] = slf2.slf.tags['times'][times2]
            for time in range(len(slf2.slf.tags['cores']))[it2:]:
                slf2.slf.tags['cores'].remove(slf2.slf.tags['cores'][-1])
            for time in range(len(slf2.slf.tags['cores']))[::-1]:
                if time not in times2:
                    slf2.slf.tags['cores'].remove(slf2.slf.tags['cores'][time])
            times0 = list(range(len(slf2.slf.tags['times'])))
            if options.head:
                print(mes)
        else:
            comparable = False
            print('\n  /!\\ no common time frames. '\
                    'The files are not comparables.\n')
        times0 = list(range(len(slf1.slf.tags['times']))) # ... for instance
    if found and options.head:
        print(mes)

    # ~~> File geometries
    mes = ''
    if slf1.slf.npoin2 != slf2.slf.npoin2:
        mes = mes + '     + npoin2 = '+str(slf1.slf.npoin2)+' in '+slf_file1
        mes = mes + '     * npoin2 = '+str(slf2.slf.npoin2)+' in '+slf_file2
        mes = mes + '\n'
    if slf1.slf.nplan != slf2.slf.nplan:
        mes = mes + '     + nplan = '+str(slf1.slf.nplan)+' in '+slf_file1
        mes = mes + '     * nplan = '+str(slf2.slf.nplan)+' in '+slf_file2
        mes = mes + '\n'
    if mes != '':
        if options.head:
            print('\n  <> Gemetry:\n'+mes)
        comparable = False
        print('\n  /!\\different geometries. The files are not comparables.\n')

    if options.head:
        # ~~> File trangulations
        diff = slf1.slf.ikle2 - slf2.slf.ikle2
        if np.argwhere(diff > [0, 0, 0]):
            print('\n  <> 2D Triangulation:\n')
            print('     + number of mismatches: '+\
                    repr(len(np.argwhere(diff == [0, 0, 0]).T[0])))
            print('     + mismatched elements: '+\
                    repr(np.argwhere(diff == [0, 0, 0]).T[0][::3]))

    if options.head:
        # ~~> File geo-localisation
        diff = np.sqrt(np.power((slf1.slf.meshx-slf2.slf.meshx), 2) + \
                        np.power((slf1.slf.meshy-slf2.slf.meshy), 2))
        points = np.argwhere(diff > float(options.epsilon)).ravel()
        if points:
            print('\n  <> Geo-Localisation:\n')
            print('     + maximum distance between points : '+\
                  str(max(diff[points])))
            pt_th = 100*len(points)/len(diff)
            print('     + number of points above clocelyness threshold : '+\
                  str(len(points))+' (i.points. '+str(pt_th)+'% of points)')
            print('     + node numbers : '+
                  repr(np.arange(slf1.slf.npoin3)[points]))

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Core differences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if options.core and comparable:
        print('\n\nCore differences: \n'+72*'~'+'\n')

        found = False
        for time in times0:
            for var in cmn_vars:
                ivar = slf1.slf.varnames.index(var)
                jvar = slf2.slf.varnames.index(var)
                var1 = slf1.slf.get_variables_at(time,
                                                 [slf1.slf.varindex[ivar]])
                var2 = slf2.slf.get_variables_at(time,
                                                 [slf2.slf.varindex[jvar]])
                diff = np.absolute(var1 - var2).ravel()
                points = np.argwhere(diff > float(options.epsilon)).ravel()
                if points.size != 0:
                    found = True
                    time1 = slf1.slf.tags['times'][time]
                    time2 = slf2.slf.tags['times'][time]
                    print('\n  <> Frame: '+str(time)+' (times: '+\
                          '{0:.2f}'.format(time1)+' / '+'{0:.2f}'.format(time2)+
                          '), Variable: '+var+'\n')
                    print('     + max difference: ', max(diff[points]))
                    print('     + number of values above threshold : '+
                          str(len(points))+' (i.points. '+
                          str(100*len(points)/len(diff))+
                          '% of points)')
                    print('     + node numbers :          '+
                          repr(np.arange(slf1.slf.npoin3)[points]))
                    print('     + values at those nodes : '+
                          repr(diff[np.arange(slf1.slf.npoin3)[points]]))
        if not found:
            print('  <> None to the epsilon: '+repr(options.epsilon))

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)

if __name__ == "__main__":
    main()
