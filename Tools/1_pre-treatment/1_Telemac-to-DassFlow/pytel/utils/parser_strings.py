r"""
   @note ... this work is based on a collaborative effort
    of the Telemac-Mascaret consortium

   @history 15/11/2011 -- Sebastien E. Bourban

   @history 25/08/2013 -- Sebastien E. Bourban and Juliette C. Parisi
      Complete re-work of the definition of points and frames:
      - points can include 3D points with vairous vertical references to planes
      - frames can include ranges
      These are mainly used to parse the keys "extract" and "time" set in the
         XML files for the validation or the extraction of data.

   @history 23/09/2014 -- Sebastien E. Bourban
      parse_array_grid has been split to include both 2D and 3D grids

   @brief
         Various method to parse strings into values, arrays, etc.
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import re
import sys
try:
    # Available in python 3.5
    from math import gcd
except ImportError:
    from fractions import gcd
from utils.exceptions import TelemacException

# _____                  ___________________________________________
# ____/ General Toolbox /__________________________________________/
#

RNG2D = re.compile(r'(?P<n>[\d:+-]+)', re.I)
NOD2D = re.compile(r'(?P<n>\d+)', re.I)
PNT2D = re.compile(r'\((?P<n>[^\(\)]+?)\)', re.I)
SPL2D = re.compile(r'\{(?P<n>[^\(\)]+?)\}', re.I)

EMPTY = re.compile(r'[\(\[][\)\]]', re.I)
PLANS = re.compile(r'\[(?P<n>[\d;,]+?)\]', re.I)
NUMBR = re.compile(r'(?P<n>[\d.+-dDeE]+?)')

SIMPLE = re.compile(r'(?P<n>([\d:+-]+|\([\d.+-dDeE]+?\)))')
COMPLX = re.compile(
    r'(?P<n>(\d+|\{[\d.+-dDeE]+?[;,][\d.+-dDeE]+?\}|\([\d.+-dDeE]+?[;,][\d.+-dDeE]+?\)|[\[\(][\]\)])((?=[^#@])|[#@][\[\(].+?[\]\)]|[#@][^,;\[\(\)\]]+))',
    re.I)

SQUOTE = re.compile(r"(?P<squot>'.*?')")  # ,re.I)
DQUOTE = re.compile(r'(?P<dquot>".*?")')  # ,re.I)


# _____                       ______________________________________
# ____/ General Time Jumping /_____________________________________/
#

def parse_array_frame(string, size=-1):
    """
    @brief     Decoding structure all in order
       The list of frames is delimiting points either by ',' or ';',
          and the ranges by ':'
       The output is an arry [..]. Each term is either:
          (a) an integer, representing a frame or a node or a plane for instance
          (b) a 1D-tuple of a real value, representing a time or a depth
          (c) a 3D-tuple of integers, representing an array range [0;-1;1] by
              default

    @param string (str) string to parse
    @param size (int) size

    examples of input / output
       '5'         =>  [5]
       '[4]'       =>  [4]
       '[5,6,7,0]' =>  [5, 6, 7, 0]
       '(5.6)'     =>  [(5.6,)]
       '(76);(4),[(3.3);4:14:2;0:6;8]'
                   =>  [(76.0,), (4.0,), (3.3,), (4, 14, 2), (0, 6, 1), 8]

    """
    frames = []

    # ~~ Special deal of all times ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if string == '[]':
        if size >= 0:
            return [range(size)]
        else:
            return [[0, -1, 1]]

    # ~~ Identify individual frames ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for node in re.findall(SIMPLE, string + ','):

        # ~~> Is it a time (itime) or a frame / range
        itime = node[1]
        proci = re.match(RNG2D, itime)
        procr = re.match(PNT2D, itime)
        if proci:
            r_t = proci.group('n').split(':')
            if len(r_t) == 1:
                frame_a = [int(r_t[0])]
                if size >= 0:
                    if frame_a[0] < 0:
                        frame_a[0] = max(0, size + frame_a[0])
                    else:
                        frame_a[0] = min(frame_a[0], size - 1)
                    frame_a = range(frame_a[0], frame_a[0] + 1, 1)
            else:
                if len(r_t) == 2:
                    frame_a = [int(r_t[0]), int(r_t[1]), 1]
                if len(r_t) == 3:
                    frame_a = [int(r_t[0]), int(r_t[1]), int(r_t[2])]
                if size >= 0:
                    if frame_a[0] < 0:
                        frame_a[0] = max(0, size + frame_a[0])
                    else:
                        frame_a[0] = min(frame_a[0], size - 1)
                    if frame_a[1] < 0:
                        frame_a[1] = max(0, size + frame_a[1])
                    else:
                        frame_a[1] = min(frame_a[1], size - 1)
                    frame_a = range(frame_a[0], frame_a[1] + 1, frame_a[2])
        elif procr:
            frame_a = (float(procr.group('n')),)
        else:
            raise TelemacException(\
                    '... could not parse the point <{}>'
                    ' from the string "{}"'.format(node[0], string))

        # ~~> Final packing
        frames.extend(frame_a)

    return frames


# _____                        _____________________________________
# ____/ General Space Jumping /____________________________________/
#

def parse_array_point(string, size=-1):
    """
    @brief     Decoding structure all in order
       The list of frames is delimiting points either by ',' or ';',
          and the ranges by ':'
       The output is an arry [..]. Each term is complicated ...

    examples of input / output
       # either a 2D node value or a vertical 1D profile covering all planes
       # above the 2D node
       '5'  =>  [(5, [(0, -1, 1)])]
       '(5)' =>  [(5, [(0, -1, 1)])]
       '9@2,58#3,18,4#1,4#1,76@0.e-3,8@0.5'
          =>  [(9, ([2.0, -1],)),
               (58, [3]),
               (18, [(0, -1, 1)]),
               (4, [1]),
               (4, [1]),
               (76, ([0.0, -1],)),
               (8, ([0.5, -1],))]
       '(4,5,6),[]#900'
          =>  [((4.0, 5.0, 6.0), [(0, -1, 1)]),
               ([], [900])]
       '(3;4,5)#[]'
          =>  [(3, [(0, -1, 1)]),
               (4, [(0, -1, 1)]),
               (5, [(0, -1, 1)])
       '(4;5,6)#[5:4;6;0:-1:2]'
          =>  [((4.0, 5.0, 6.0), [(5, 4, 1), 6, (0, -1, 2)])]
       '9@2,58#3,18,(4;7)#1,4#1,(76;4)@1.e-1,[8@(0.5;0.7)'
          =>  [(9, ([2.0, -1],)),
               (58, [3]),
               (18, [(0, -1, 1)]),
               ((4.0, 7.0), [1]),
               (4, [1]),
               ((76.0, 4.0), ([0.1, -1],)),
               (8, ([0.5, -1],[0.7, -1]))]
       '(4;5,6)#[5;6]'
          =>  [((4.0, 5.0, 6.0), [5, 6])]
       '(4;5,6)@(-5#3;6)'
          =>  [((4.0, 5.0, 6.0), ([-5.0, 3], [6.0, -1]))]
    """

    points = []

    # ~~ Special deal of all points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if string == '':
        if size >= 0:
            return [([], range(size))]
        else:
            return [([], [0])]

    # ~~ Identify individual points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for node in re.findall(COMPLX, string + ','):

        # ~~> Is it a point (x,y) or a node n
        x_y = node[1]
        proci = re.match(NOD2D, x_y)
        procs = re.match(SPL2D, x_y)
        procr = re.match(PNT2D, x_y)
        proce = re.match(EMPTY, x_y)
        if proci:
            point_a = int(proci.group('n'))
        elif procr:
            x_y = procr.group('n').replace(',', ';').split(';')
            if len(x_y) == 2:
                point_a = (float(x_y[0]), float(x_y[1]))
            if len(x_y) != 2:
                raise TelemacException(\
                        '... we are not allowing anything '
                        'but a pair (x,y). You have: <{}> '
                        'from the string "{}"\n'
                        '+> if you need (x,y,z) you should use a depth '
                        'above plan 0: (x,y)@z#0'.format(node[0], string))
        elif proce:
            point_a = []
        elif procs:
            x_y = procs.group('n').replace(',', ';').split(';')
            if len(x_y) == 2:
                point_a = (int(x_y[0]), int(x_y[1]))
            elif len(x_y) == 3:
                point_a = (int(x_y[0]), int(x_y[1]), int(x_y[2]))
            else:
                raise TelemacException(\
                  '... could not parse the number of re-sampling steps. '
                  'You have: <{}>from the string "{}"'.format(node[0], string))

            points.append(point_a)
            continue
        else:
            raise TelemacException(\
                    '... could not parse the point <{}> '
                    'from the string "{}"'.format(node[0], string))

        # ~~> Is it a depth d or a plane p or both
        point_b = []
        if node[2] != '':
            t_p = node[2][0]
            z_p = node[2][1:]
            if t_p == '#':  # this is a plane or a series of planes
                proci = re.match(RNG2D, z_p)
                if proci:
                    z_p = '[' + z_p + ']'
                point_b = parse_array_frame(z_p, size)
            # this is a depth or a series of depth, referenced by planes
            if t_p == '@':
                procr = re.match(NUMBR, z_p)
                if procr:
                    z_p = '(' + z_p + ')'
                procp = re.match(PNT2D, z_p)
                if procp:
                    point_b = []
                    for pnt in procp.group('n').replace(',', ';').split(';'):
                        if '#' in pnt:
                            tmp_a, tmp_b = pnt.split('#')
                            point_b.append([float(tmp_a), int(tmp_b)])
                        else:
                            # from the surface plane by default
                            point_b.append([float(pnt), -1])
                    point_b = tuple(point_b)
        else:
            if size >= 0:
                point_b = range(size)
            else:
                point_b = [0, -1, 1]

        # ~~> Final packing
        points.append((point_a, point_b))

    return points


def parse_array_grid(string, size):
    """
    @brief     Decoding structure all in order
       The grid is defined by two points and an array of re-sampling steps
       The input 'size' is either:
          - in 2D a pair of 2D points ( bottom-left, top-right )
          - in 3D a pair of 2D points and a range of planes
       The input 'size' is a pair of complex points (2D or 3D) and
          a set of re-sampling numbers
       The output is an arry [..]. Each term is complicated ...
    """

    grids = []
    minz = 0.
    maxz = 0.
    minp = 0
    maxp = 0

    if len(size) == 3:
        (minx, miny), (maxx, maxy), (minp, maxp) = size
    elif len(size) == 2:
        if len(size[0]) == 2:
            (minx, miny), (maxx, maxy) = size
        else:
            (minx, miny, minz), (maxx, maxy, maxz) = size
    n_z = maxp - minp

    # ~~ Special deal of all points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if string == '[]':
        # arbitrary value of 20 points
        dist = (maxy - miny + maxx - minx) / 20.0
        dist = min(dist, maxx - minx)
        dist = min(dist, maxy - miny)
        x_o = (maxx + minx) / 2.0
        y_o = (maxy + miny) / 2.0
        n_x = max(2, int((maxx - minx) / dist))
        n_y = max(2, int((maxy - miny) / dist))
        dist = min(dist, (maxx - minx) / (1.0 * n_x))
        dist = min(dist, (maxy - miny) / (1.0 * n_y))
        if len(size) == 2 and len(size[0]) == 2:
            return [[(x_o - n_x * dist / 2.0, y_o - n_y * dist / 2.0),
                     (x_o + n_x * dist / 2.0, y_o + n_y * dist / 2.0),
                     [n_x, n_y]]]
        # TODO: make sure you can suport this option
        elif len(size) == 2 and len(size[0]) == 3:
            z_o = (maxz + minz) / 2.0
            n_z = 10
            dizt = (maxx - minx) / (1.0 * n_z)  # arbitrary value of 10 points
            return [[(x_o - n_x * dist / 2.0, y_o - n_y * dist / 2.0, z_o \
                      - n_z * dizt / 2.0),
                     (x_o + n_x * dist / 2.0, y_o + n_y * dist / 2.0, z_o\
                      + n_z * dizt / 2.0),
                     [n_x, n_y, n_z]]]
        else:
            return [[(x_o - n_x * dist / 2.0, y_o - n_y * dist / 2.0),
                     (x_o + n_x * dist / 2.0, y_o + n_y * dist / 2.0),
                     range(minp, maxp),
                     [n_x, n_y, n_z]]]

    # ~~ Decoding of user entrance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    s2g = parse_array_point(string)
    if gcd(len(s2g), 3) != 3:
        raise TelemacException(\
                '... could not parse your grid . "{}."'
                'It should be triplets made of 2 points (;)(;) '
                'and an array of resampling steps {{;}}.'.format(string))

    for i in range(len(s2g) // 3):
        pta, ptb, n_p = s2g[3 * i:3 * (i + 1)]
        if len(n_p) == 2:
            grids.append([pta[0], ptb[0], n_p])
        # TODO: support a range of fixed depths as well as fixed planes
        elif len(n_p) == 3:
            z_p = '[' + str(pta[1][0]) + ':' + str(ptb[1][0]) + ']'
            grids.append([pta[0], ptb[0], parse_array_frame(z_p, n_z), n_p])

    return grids


def parse_array_paires(string):
    paires = []

    # ~~ Special deal of all paires ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if string == '':
        return []

    # ~~ Identify individual paires ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for node in re.findall(COMPLX, string + ','):

        # ~~> Is it a point (x,y) or a node n
        x_y = node[1]
        proci = re.match(NOD2D, x_y)
        procr = re.match(PNT2D, x_y)
        proce = re.match(EMPTY, x_y)
        if proci:
            point_a = int(proci.group('n'))
        elif procr:
            x_y = procr.group('n').replace(',', ';').split(';')
            if len(x_y) == 2:
                point_a = (float(x_y[0]), float(x_y[1]))
            if len(x_y) != 2:
                raise TelemacException(\
                    '... we are not allowing anything anything '
                    'but a pair (x,y). You have: <{}> from the string "{}"\n'
                    '   +> if you need (x,y,z) you should use a depth above '
                    'plan 0: (x,y)@z#0'.format(node[0], string))
        elif proce:
            point_a = []
        else:

            raise TelemacException(\
                    '... could not parse the point <{}>'
                    'from the string "{}"'.format(node[0], string))
        # ~~> Final packing
        paires.append(point_a)

    return paires


# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

__author__ = "Sebastien E. Bourban"
__date__ = "$15-Nov-2011 08:51:29$"


def main():
    # ~~ space ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if True:  # S.E.Bourban completed testing on August 25, 2013
        print('\n\n')
        strings = ['5', '9@2,58#3,18,4#1,4#1,76@0.e-3,8@0.5', '(4,5,6),[]#900',
                   '(3;4,5)#[]', '(4;5,6)#[5:4;6;0:-1:2]',
                   '(5)', '9@2,58#3,18,(4;7)#1,4#1,(76;4)@1.e-1,[8@(0.5;0.7)',
                   '(4;5,6)#[5;6]', '(4;5,6)@(-5#3;6)']
        for string in strings:
            print(string + ' => ' + repr(parse_array_point(string)))

    # ~~ time ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if True:  # S.E.Bourban completed testing on August 25, 2013
        print('\n\n')
        strings = ['5', '[4]', '[5,6,7,0]',
                   '(5.6)', '(76);(4),[(3.3);4:14:2;0:6;8]']
        for string in strings:
            print(string + ' => ' + repr(parse_array_frame(string)))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)


if __name__ == "__main__":
    main()
