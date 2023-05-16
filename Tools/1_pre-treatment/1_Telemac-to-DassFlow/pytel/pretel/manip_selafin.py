r"""@author Sebastien E. Bourban

"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
from os import path
import numpy as np
# ~~> dependencies towards other modules

from compilation.parser_fortran import clean_quotes
from data_manip.extraction.parser_lqd import LQD
from data_manip.extraction.parser_kenue import InS
from data_manip.conversion import convert_utm as utm
from data_manip.formats.selafins import Selafins
from data_manip.formats.selafin import Selafin
from utils.files import move_file
from utils.exceptions import TelemacException

from pretel.meshes import tessellate_poly
from pretel.scan_selafin import ScanSelafin
from pretel.alter_selafin import AlterSelafin
from pretel.chop_selafin import ChopSelafin
from pretel.crunch_selafin import CrunchSelafin
from pretel.scan_spectral import ScanSpectral
from pretel.sub_selafin import SubSelafin
from pretel.calcs_selafin import CalcsSelafin
from pretel.transf_selafin import TransfSelafin

def scan(options):
    """
    Scan of a file
    """
    slf_files = options.input_files
    for slf_file in slf_files:
        slf_file = path.realpath(slf_file)
        if not path.exists(slf_file):
            raise TelemacException(\
                    '\nCould not find the file named: {}'.format(slf_file))
        print('\n\nScanning ' + path.basename(slf_file) + ' within ' + \
                path.dirname(slf_file) + '\n'+'~'*72+'\n')
        vrs = options.xvars
        if options.xvars is not None:
            vrs = clean_quotes(options.xvars.replace('_', ' '))
        times = (int(options.tfrom), int(options.tstep), int(options.tstop))
        slf = ScanSelafin(slf_file, times=times, vrs=vrs)
        slf.print_header()
        if options.core:
            slf.print_core()
        else:
            slf.print_time_summary()


def spec(options):
    """
    Spectral file
    """
    slf_files = options.input_files
    for slf_file in slf_files:

        slf_file = path.realpath(slf_file)
        if not path.exists(slf_file):
            raise TelemacException(\
                 '\nCould not find the file named: {}'.format(slf_file))
        print('\n\nScanning ' + path.basename(slf_file) + ' within ' + \
                path.dirname(slf_file) + '\n'+'~'*72+'\n')
        vrs = options.xvars
        if options.xvars is not None:
            vrs = clean_quotes(options.xvars.replace('_', ' '))
        times = (int(options.tfrom), int(options.tstep), int(options.tstop))
        slf = ScanSpectral(slf_file, times=times, vrs=vrs)
        slf.print_header()
        if options.core:
            slf.print_core(int(options.accuracy))
        else:
            slf.print_time_summary()


def chop(options):
    """
    Chopping of a file
    """
    root_file = None
    if not options.freplace:
        if not options.parallel:
            if len(options.args) != 2:
                raise TelemacException(\
                        '\nThe code "chop" (without --replace) '
                        'here requires 2 file names\n')
            slf_files = [options.args[0]]
            out_file = options.args[1]
        else:
            if len(options.args) != 3:
                raise TelemacException(\
                        '\nThe code "chop" (without --replace) '
                        'here requires 2 file names and '
                        '1 file root name for the partition\n')
            slf_files = [options.args[0]]
            root_file = options.args[1]
            out_file = options.args[2]
    else:
        slf_files = options.args
        out_file = "chop-tmp.slf"

    for slf_file in slf_files:

        slf_file = path.realpath(slf_file)
        if not path.exists(slf_file):
            raise TelemacException(\
               '\nCould not find the file named: {}'.format(slf_file))
        print('\n\nChoping ' + path.basename(slf_file) + ' within ' + \
                path.dirname(slf_file) + '\n'+'~'*72+'\n')
        vrs = options.xvars
        if options.xvars is not None:
            vrs = clean_quotes(options.xvars.replace('_', ' '))
        times = (int(options.tfrom), int(options.tstep), int(options.tstop))
        slf = ChopSelafin(slf_file, times=times, vrs=vrs, root=root_file)
        if options.eswitch:
            slf.alter_endian()
        if options.fswitch:
            slf.alter_float()

        slf.put_content(out_file)

        if options.freplace:
            move_file(out_file, slf_file)


def alter(options):
    """
    Modifications in the file
    """
    root_file = None
    if not options.freplace:
        if not options.parallel:
            if len(options.args) != 2:
                raise TelemacException(\
                        '\nThe code "alter" (without --replace) '
                        'requires 2 file names\n')
            slf_files = [options.args[0]]
            out_file = options.args[1]
        else:
            if len(options.args) != 3:
                raise TelemacException(\
                        '\nThe code "alter" (without --replace) '
                        'here requires 2 file names and '
                        '1 file root name for the partition\n')
            slf_files = [options.args[0]]
            root_file = options.args[1]
            out_file = options.args[2]
    else:
        slf_files = options.args
        out_file = "chop-tmp.slf"

    for slf_file in slf_files:
        slf_file = path.realpath(slf_file)
        if not path.exists(slf_file):
            raise TelemacException(\
                '\nCould not find the file named: {}'.format(slf_file))
        print('\n\nAltering ' + path.basename(slf_file) + ' within ' + \
                path.dirname(slf_file) + '\n'+'~'*72+'\n')
        vrs = options.xvars
        if options.xvars is not None:
            vrs = clean_quotes(options.xvars.replace('_', ' '))
        times = (int(options.tfrom), int(options.tstep), int(options.tstop))
        slf = AlterSelafin(slf_file, times=times, vrs=vrs, root=root_file)
        if options.atitle is not None:
            slf.alter_title(options.atitle)
        if options.areset:
            slf.alter_times(p_t=-slf.slf.tags['times'][0])
        if options.adate is not None:
            slf.alter_datetime(date=options.adate.split('-'))
        if options.atime is not None:
            slf.alter_datetime(time=options.atime.split(':'))
        if options.aswitch:
            slf.switch_vars()
        if options.eswitch:
            slf.alter_endian()
        if options.fswitch:
            slf.alter_float()
        if options.aname is not None:
            slf.alter_vars(options.aname)
        slf.alter_times(m_t=float(options.atm), p_t=float(options.atp))
        slf.alter_mesh(m_x=float(options.axm), p_x=float(options.axp),
                       m_y=float(options.aym), p_y=float(options.ayp))
        if options.azname is not None:
            slf.alter_values(options.azname,
                             m_z=float(options.azm), p_z=float(options.azp))
        if options.sph2ll is not None:
            radius = 6371000.
            long0, lat0 = options.sph2ll.split(":")
            long0 = np.deg2rad(float(long0))
            lat0 = np.deg2rad(float(lat0))
            const = np.tan(lat0/2. + np.pi/4.)
            slf.slf.meshx = np.rad2deg(slf.slf.meshx/radius + long0)
            expo = np.exp(slf.slf.meshy/radius)
            slf.slf.meshy = np.rad2deg(2.*np.arctan(const*expo) - np.pi/2.)
        if options.ll2sph is not None:
            radius = 6371000.
            long0, lat0 = options.ll2sph.split(":")
            long0 = np.deg2rad(float(long0))
            lat0 = np.deg2rad(float(lat0))
            slf.slf.meshx = radius * (np.deg2rad(slf.slf.meshx) - long0)
            slf.slf.meshy = radius * \
                      (np.log(np.tan(np.deg2rad(slf.slf.meshy)/2. + np.pi/4.)) \
                              - np.log(np.tan(lat0/2. + np.pi/4.)))
        if options.ll2utm is not None:
            if options.ll2utm != 'XXX':
                zone = int(options.ll2utm[:-1])
                zone_letter = options.ll2utm[-1]
            else:
                zone = None
                zone_letter = None
            slf.slf.meshx, slf.slf.meshy, zone, zone_letter = \
                      utm.from_latlon(slf.slf.meshx, slf.slf.meshy,
                                      force_zone_number=zone,
                                      force_zone_letter=zone_letter)
            print(' ~> Converterted to UTM ZONE {}{}'.format(zone, zone_letter))
        if options.utm2ll is not None:
            zone = int(options.utm2ll[:-1])
            zone_letter = options.utm2ll[-1]
            slf.slf.meshx, slf.slf.meshy = \
                      utm.to_latlon(slf.slf.meshx, slf.slf.meshy,
                                    zone, zone_letter)

        slf.put_content(out_file)

        if options.freplace:
            move_file(out_file, slf_file)


def merge(options):
    """
    Merging two selafin files
    """
    root_file = None
    if not options.parallel:
        if len(options.args) < 3:
            raise TelemacException(\
                    '\nThe code "merge" requires '
                    'at leat 2 file names, aside '
                    'from the options\n')
        slf_files = options.args[0:len(options.args)-1]
        out_file = options.args[len(options.args)-1]

        slfs = Selafins()
        print('\n\nMerging into ' + path.basename(out_file) + ' within ' + \
                path.dirname(out_file) + '\n'+'~'*72+'\n')
        for slf_file in slf_files:
            slf_file = path.realpath(slf_file)
            if not path.exists(slf_file):
                raise TelemacException(\
                        '\nCould not find '
                        'the file named: {}'.format(slf_file))
            slfs.add(slf_file)

        slfs.put_content(out_file)

    else:
        if len(options.args) != 3:
            raise TelemacException(\
                    '\nThe code "merge" here requires '
                    '2 file names and '
                    '1 file root name for the partition\n')
        slf_file = options.args[0]
        root_file = options.args[1]
        out_file = options.args[2]

        print('\n\nMerging into ' + path.basename(out_file) + ' within ' \
                + path.dirname(out_file) + '\n'+'~'*72+'\n')
        slf_file = path.realpath(slf_file)
        if not path.exists(slf_file):
            raise TelemacException(\
                    '\nCould not find '
                    'the file named: {}'.format(slf_file))

        vrs = options.xvars
        if options.xvars is not None:
            vrs = clean_quotes(options.xvars.replace('_', ' '))
        times = (int(options.tfrom), int(options.tstep), int(options.tstop))
        slf = AlterSelafin(slf_file, times=times, vrs=vrs, root=root_file)
        if options.atitle is not None:
            slf.alter_title(options.atitle)
        if options.areset:
            slf.alter_times(p_t=-slf.slf.tags['times'][0])
        if options.adate is not None:
            slf.alter_datetime(date=options.adate.split('-'))
        if options.atime is not None:
            slf.alter_datetime(time=options.atime.split(':'))
        if options.aswitch:
            slf.switch_vars()
        if options.eswitch:
            slf.alter_endian()
        if options.fswitch:
            slf.alter_float()
        if options.aname is not None:
            slf.alter_vars(options.aname)
        slf.alter_times(m_t=float(options.atm), p_t=float(options.atp))
        slf.alter_mesh(m_x=float(options.axm), p_x=float(options.axp),
                       m_y=float(options.aym), p_y=float(options.ayp))
        if options.azname is not None:
            slf.alter_values(options.azname,
                             m_z=float(options.azm), p_z=float(options.azp))

        slf.put_content(out_file)

def diff(options):
    """
    diff between two serafin files
    """
    if len(options.args) < 2:
        raise TelemacException(\
                '\nThe code "diff" uses a minimum of '
                '3 argumensts, aside from the options\n')
    slf_files = options.args[0:len(options.args)-1]
    out_file = options.args[len(options.args)-1]

    slfs = Selafins()
    print('\n\nDifferences into {}\n{}\n'.format(path.basename(out_file),
                                                 '~'*72))
    for slf_file in slf_files:
        slf_file = path.realpath(slf_file)
        if not path.exists(slf_file):
            raise TelemacException(\
                    '\nCould not find '
                    'the file named: {}'.format(slf_file))
        slfs.add(slf_file)
    slfs.put_content(out_file)


def sample(options):
    """
    Set liquid boundary file from a selafin file
    """
    root_file = None
    if not options.parallel:
        if len(options.args) < 3:
            raise TelemacException(\
                    '\nThe code "sample" requires '
                    'at least 2 file names and '
                    'one series of node numbers\n')
        slf_file = options.args[0]
        out_file = options.args[1]
        nod_list = []
        for nod in options.args[2].split(" "):
            nod_list.append(int(nod))
    else:
        if len(options.args) != 4:
            raise TelemacException(\
                    '\nThe code "sample" here '
                    'requires 2 file names, '
                    '1 file root name for the partition and '
                    '1 series of node numbers\n')
        slf_file = options.args[0]
        root_file = options.args[1]
        out_file = options.args[2]
        nod_list = []
        for nod in options.args[3].split(" "):
            nod_list.append(int(nod))

    slf_file = path.realpath(slf_file)
    if not path.exists(slf_file):
        raise TelemacException(\
                '\nCould not find the file named: '
                '{}'.format(slf_file))
    print('\n\nSample ' + path.basename(slf_file) + ' within ' + \
            path.dirname(slf_file) + '\n'+'~'*72+'\n')
    vrs = options.xvars
    if options.xvars is not None:
        vrs = clean_quotes(options.xvars.replace('_', ' '))
    times = (int(options.tfrom), int(options.tstep), int(options.tstop))
    slf = ChopSelafin(slf_file, times=times, vrs=vrs, root=root_file)

    lqd = LQD(vrs=[zip(slf.slf.varnames, slf.slf.varunits), nod_list],
              date=slf.slf.datetime, times=slf.slf.tags['times'],
              series=slf.get_series(nod_list))
    lqd.put_content(out_file)

def subdivide(options):
    """
    Subdivide a mesh
    """
    if not options.freplace:
        if len(options.args) != 2:
            raise TelemacException(\
                    '\nThe code "subdivide" '
                    '(without --replace) here '
                    'requires 2 file names\n')
        slf_file = options.args[0]
        out_file = options.args[1]
    else:
        if len(options.args) != 1:
            raise TelemacException(\
                    '\nThe code "subdivide" (with --replace) '
                    'here requires 1 file name at a time\n')
        slf_file = options.args[0]
        out_file = "subdivide-tmp.slf"

    slf_file = path.realpath(slf_file)
    if not path.exists(slf_file):
        raise TelemacException(\
                '\nCould not find'
                ' the file named: {}'.format(slf_file))
    print('\n\nSubdividing ' + path.basename(slf_file) + ' within ' + \
            path.dirname(slf_file) + '\n'+'~'*72+'\n')
    slf = SubSelafin(slf_file)
    slf.put_content(out_file)

    if options.freplace:
        move_file(out_file, slf_file)


def tesselate(options):
    """
    Generate a mesh from a polygon
    """
    if not options.freplace:
        if len(options.args) != 2:
            raise TelemacException(\
                    '\nThe code "tessellate" here '
                    'requires one i2s/i3s file and '
                    'one output slf file\n')
        i3s_file = options.args[0]
        out_file = options.args[1]
    else:
        if len(options.args) != 1:
            raise TelemacException(\
                    '\nThe code "tessellate" here '
                    'requires one i2s/i3s file\n')
        i3s_file = options.args[0]
        head, _ = path.splitext(i3s_file)
        out_file = head+'.slf'

    i3s_file = path.realpath(i3s_file)
    if not path.exists(i3s_file):
        raise TelemacException(\
                '\nCould not find '
                'the file named: {}'.format(i3s_file))

    print('\n\nTessellating ' + path.basename(i3s_file) + ' within ' + \
            path.dirname(i3s_file) + '\n'+'~'*72+'\n')
    i2s = InS(i3s_file)
    ikle2, ipob2, meshx, meshy = tessellate_poly(i2s, debug=True)

    print('\n\nWriting down the Selafin file ' + \
            path.basename(out_file) + '\n'+'~'*72+'\n')
    slf = Selafin('')
    slf.title = ''
    slf.nplan = 1
    slf.ndp2 = 3
    slf.ndp3 = 3
    slf.nbv1 = 1
    slf.nvar = 1
    slf.varindex = 1
    slf.varnames = ['BOTTOM          ']
    slf.varunits = ['M               ']
    slf.ikle2 = ikle2
    slf.ikle3 = slf.ikle2
    slf.meshx = meshx
    slf.meshy = meshy
    slf.npoin2 = i2s.npoin
    slf.npoin3 = slf.npoin2
    slf.nelem2 = len(slf.ikle2)/slf.ndp3
    slf.nelem3 = slf.nelem2
    slf.iparam = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
    slf.ipob2 = ipob2
    slf.ipob3 = slf.ipob2
    slf.fole = {'hook':open(out_file, 'wb'), 'endian':">",
                'float':('f', 4), 'name':out_file}
    slf.tags['times'] = [1]
    if options.sph2ll is not None:
        radius = 6371000.
        long0, lat0 = options.sph2ll.split(":")
        long0 = np.deg2rad(float(long0))
        lat0 = np.deg2rad(float(lat0))
        const = np.tan(lat0/2. + np.pi/4.)
        slf.meshx = np.rad2deg(slf.meshx/radius + long0)
        slf.meshy = np.rad2deg(2.*np.arctan(const*np.exp(slf.meshy/radius)) \
                                      - np.pi/2.)
    if options.ll2sph is not None:
        radius = 6371000.
        long0, lat0 = options.ll2sph.split(":")
        long0 = np.deg2rad(float(long0))
        lat0 = np.deg2rad(float(lat0))
        slf.meshx = radius * (np.deg2rad(slf.meshx) - long0)
        slf.meshy = radius * \
                  (np.log(np.tan(np.deg2rad(slf.meshy)/2. + np.pi/4.)) \
                                      - np.log(np.tan(lat0/2. + np.pi/4.)))
    if options.ll2utm is not None:
        if options.ll2utm != 'XXX':
            zone = int(options.ll2utm[:-1])
            zone_letter = options.ll2utm[-1]
        else:
            zone = None
            zone_letter = None
        slf.meshx, slf.meshy, zone, zone_letter = \
                  utm.from_latlon(slf.meshx, slf.meshy,
                                  force_zone_number=zone,
                                  force_zone_letter=zone_letter)
        print(' ~> Converterted to UTM ZONE {}{}'.format(zone, zone_letter))
    if options.utm2ll is not None:
        zone = int(options.utm2ll[:-1])
        zone_letter = options.utm2ll[-1]
        slf.meshx, slf.meshy = \
                  utm.to_latlon(slf.meshx, slf.meshy,
                                zone, zone_letter)
    slf.append_header_slf()
    slf.append_core_time_slf(0)
    slf.append_core_vars_slf([np.zeros(slf.npoin2)])
    slf.fole['hook'].close()

def calcs(options, code_name):
    """
    Doing calcs, crunh, transf
    """
    root_file = None
    if not options.parallel:
        if len(options.args) < 2:
            raise TelemacException(\
                    '\nThe code "calcs" requires 2 file names\n')
        slf_file = options.args[0]
        out_file = options.args[1]
    else:
        if len(options.args) != 3:
            raise TelemacException(\
                    '\nThe code "calcs" requires '
                    '2 file names and 1 root file name '
                    'for parallel inputs\n')
        slf_file = options.args[0]
        root_file = options.args[1]
        out_file = options.args[2]

    slf_file = path.realpath(slf_file)
    if not path.exists(slf_file):
        raise TelemacException(\
                '\nCould not find the file named: {}'.format(slf_file))
    print('\n\nCalculations for ' + path.basename(slf_file) + ' within ' + \
            path.dirname(slf_file) + '\n'+'~'*72+'\n')
    vrs = options.xvars
    calc_list = []
    if options.xvars is not None:
        vrs = clean_quotes(options.xvars.replace('_', ' '))
        calc_list = vrs.split(':')
    if code_name == 'calcs':
        times = (int(options.tfrom), int(options.tstep), int(options.tstop))
        slf = CalcsSelafin(slf_file, times=times, root=root_file)
        print('   ~> Assembling the following variables together '\
                'into the file:')
        for calc in calc_list:
            if calc.upper() in "WATER DEPTH":
                print('      +> WATER DEPTH')
                slf.calc_water_depth()
            if calc.upper() in "KINETIC ENERGY":
                print('      +> KINETIC ENERGY')
                slf.calc_kinetic_energy()
    elif code_name == 'transf':
        times = (float(options.tfrom), float(options.tstep),
                 float(options.tstop))
        slf = TransfSelafin(slf_file, times=times, root=root_file)
        print('   ~> Computing an animation for the following variable(s):')
        for calc in calc_list:
            if calc.upper() in "WAVE SURFACE":
                print('      +> WAVE SURFACE')
                slf.calc_free_surface_from_artemis()
    elif code_name == 'crunch':
        times = (int(options.tfrom), int(options.tstep), int(options.tstop))
        slf = CrunchSelafin(slf_file, times=times, root=root_file)
        print('   ~> Assembling the following variables into the file:')
        for calc in calc_list:
            if calc.upper() in "SURFACE RANGE":
                print('      +> SURFACE RANGE')
                slf.calc_surface_range()
            if calc.upper() in "MAXIMUM SPEED":
                print('      +> MAXIMUM SPEED')
                slf.calc_maximum_speed()
            if calc.upper() in "TIME OF PEAK":
                print('      +> TIME OF PEAK')
                slf.calc_peak_time_modulo_m2()
            if calc.upper() in "RESIDUAL U":
                print('      +> RESIDUAL U')
                slf.calc_residual_velocity()

    slf.alter_times(m_t=float(options.atm), p_t=float(options.atp))
    slf.alter_mesh(m_x=float(options.axm), p_x=float(options.axp),
                   m_y=float(options.aym), p_y=float(options.ayp))
    if options.azname is not None:
        slf.alter_values(options.azname,
                         m_z=float(options.azm), p_z=float(options.azp))
    if options.eswitch:
        slf.alter_endian()
    if options.fswitch:
        slf.alter_float()

    slf.put_content(out_file)
