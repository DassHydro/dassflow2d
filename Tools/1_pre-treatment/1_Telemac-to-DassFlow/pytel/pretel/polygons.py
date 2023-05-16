r"""@author Sebastien E. Bourban

    @history 26/12/2011 -- Sebastien E. Bourban

    @brief
        Tools for trivial polygon operations
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
from os import path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import numpy as np
sys.path.append(path.join(path.dirname(sys.argv[0]), '..'))
# ~~> dependencies towards other modules
# ~~> dependencies towards other pytel/modules
from utils.geometry import get_cone_angle, is_close, get_norm2
from data_manip.extraction.parser_kenue import InS

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#

# _____                  ___________________________________________
# ____/ General Toolbox /__________________________________________/
#
def join_segments(poly_lines):

    poly_gones = []
    while poly_lines != []:
        # ~~> starting point
        poly = poly_lines[0]
        i, j = poly[0], poly[len(poly)-1]
        # ~~> case of closed line
        if i == j:
            poly_gones.append(poly[1:])
            poly_lines.pop(0)
            continue
        # ~~> iterative process
        iline = 1
        polyi = poly_lines[iline]
        while j != polyi[0]:
            iline += 1
            polyi = poly_lines[iline]
        # ~~> merging the two segments
        poly.extend(polyi[1:])
        poly_lines[0] = poly
        poly_lines.pop(iline)

    return poly_gones

def smooth_subdivise(poly, vals, p_type, weight):

    ptwice = np.zeros((2*len(poly)-1+p_type, 2))
    #vtwice = np.zeros((2*len(poly)-1+p_type,2,len(lavs)))
    # ~~> save original set
    for i in range(len(poly)):
        ptwice[2*i] = poly[i]
    # ~~> include intermediates
    for i in range(len(poly)-1):
        ptwice[2*i+1] = (poly[i]+poly[i+1])/2.
    if p_type != 0:
        ptwice[2*len(poly)-1] = (poly[0]+poly[len(poly)-1])/2.
    # ~~> weighted-average of the original
    for i in range(len(poly)-1)[1:]:
        ptwice[2*i] = weight*ptwice[2*i] + \
                          (1-weight)*(ptwice[2*i-1]+ptwice[2*i+1])/2.
    if p_type != 0:
        ptwice[0] = weight*ptwice[0] + \
                        (1-weight)*(ptwice[len(ptwice)-1]+ptwice[1])/2.
        ptwice[len(ptwice)-2] = weight*ptwice[len(ptwice)-2] + \
                                (1-weight)*(ptwice[len(ptwice)-1]+\
                                ptwice[len(ptwice)-3])/2.

    return ptwice, vals, p_type

def remove_duplicates(poly, p_type): # /!\ does not work anymore
    found = True
    while found:
        i = 0
        found = False
        while i < len(poly)-1:
            if is_close(poly[i], poly[i+1], size=10):
                found = True
                poly = np.delete(poly, i+1, 0)
            i += 1
    if len(poly) == 1:
        return [], 0
    elif len(poly) == 2:
        return [], 0 #poly,0
    else:
        if p_type != 0:
            if is_close(poly[len(poly)-1], poly[0], size=10):
                poly = np.delete(poly, len(poly)-1, 0)
        if len(poly) < 3:
            return [], 0 #poly,0
        return poly, p_type

def subsample_distance(poly, vals, p_type, dist):

    found = True
    while found:
        i = 0
        found = False
        while i < len(poly)-1:
            if dist > get_norm2(poly[i], poly[i+1]):
                poly[i] = (poly[i]+poly[i+1])/2.
                poly = np.delete(poly, i+1, 0)
                vals[i] = (vals[i]+vals[i+1])/2.
                vals = np.delete(vals, i+1, 0)
                found = True
            i += 1
    if len(poly) == 1:
        return [], [], 0
    elif len(poly) == 2:
        return [], [], 0 #poly,0
    else:
        if p_type != 0:
            if dist > get_norm2(poly[len(poly)-1], poly[0]):
                poly[len(poly)-1] = (poly[len(poly)-1]+poly[0])/2.
                poly = np.delete(poly, 0, 0)
                vals[len(vals)-1] = (vals[len(vals)-1]+vals[0])/2.
                vals = np.delete(vals, 0, 0)
        if len(poly) < 3:
            return [], [], 0 #poly,0
        return poly, vals, p_type

def subsample_angle(poly, vals, p_type, angle):

    found = True
    while found:
        i = 0
        found = False
        while i < len(poly)-4:
            if angle > 180*abs(abs(get_cone_angle(poly[i], poly[i+1],
                                                  poly[i+2]))\
                                     - np.pi)/np.pi:
                poly = np.delete(poly, i+1, 0)
                vals = vals.pop(i+1)
                found = True
            if angle > 180*abs(abs(get_cone_angle(poly[i+1], poly[i+2],
                                                  poly[i+3]))\
                               - np.pi)/np.pi:
                poly = np.delete(poly, i+2, 0)
                vals = vals.pop(i+2)
                found = True
            i += 2
    if len(poly) < 3:
        return [], [], 0 #poly,vals,0
    return poly, vals, p_type

def is_clockwise(poly):
    # assumes that poly does not duplicate points
    wise = 0
    for i in range(len(poly)):
        z = (poly[(i+1)%len(poly)][0]-poly[i][0]) \
            *(poly[(i+2)%len(poly)][1]-poly[(i+1)%len(poly)][1]) \
            - (poly[(i+1)%len(poly)][1]-poly[i][1]) \
            * (poly[(i+2)%len(poly)][0]-poly[(i+1)%len(poly)][0])
        if z > 0:
            wise += 1
        elif z < 0:
            wise -= 1
    return wise < 0

def make_clockwise(poly, vals):
    # TODO -- you need to also flip the values
    if not is_clockwise(poly):
        return np.flipud(poly), vals[::-1]
    return poly, vals

def make_anti_clockwise(poly, vals):
    # TODO -- you need to also flip the values
    if is_clockwise(poly):
        return np.flipud(poly), vals[::-1]
    return poly, vals

def get_area(poly):
    # assumes that poly does not duplicate points
    x, y = poly.T
    area = np.sum(x[:-1]*y[1:] - x[1:]*y[:-1])
    return abs(area + (x[-1]*y[0] - x[0]*y[-1]))/2.0

# _____                  ___________________________________________
# ____/ Primary Classes /__________________________________________/
#

class Polygons(object):
    """
        Polygons is a lightweight structure to hold polygon information,
            whether to support i2s/i3s or shape or other file types.
        Polygons repalces the previous object InS
        Note that the derived class of Polygons have to implement the methods
            parse_content and put_content
    """

    def __init__(self):
        self.object = None
        # important for coordinate conversions and distance calculations
        self.coordinates = {'type':None}

        # list of several polygons, each being defined as a pairs of x,y
        # self.poly = []
        # list of several polygons, each being defined as the values for the
        # corresponding nodes
        # self.vals = []
        # an tuples of integers for each polygon, where for the
        # first value: 0 = open; 1 = closed clockwise;
        #              2 = closed anti-clockwise; 3 = ...
        # second value: 0 = soft line; 1 = hard line; 2 = ...
        # self.type = []
        # a list of attributs for each polygon, common to all nodes on the
        # polygon
        # self.atrbut = []
        self.npoly = 0
        # /!\ the poly does not duplicate the first and last node for closed
        # contours
        self.npoin = 0

    def parse_content(self, file_name):
        # file parsing is based on the name of the extension
        _, tail = path.splitext(file_name)
        # ~~> Case of a Kenue type i2s/i3s file
        if tail in ['.i2s', '.i3s']:
            self.object = InS(file_name)
        self.npoly = self.object.npoly
        # ~~> Case of a Arc-GIS type shap file (which comes with additional
        # related files)
        #elif tail == '.shp':
        #   head, fileType, self.npoin, self.poly, self.vals, self.type, \
        #       self.atrbut = getShp(self.file_name)
        # ~~> Sort out the poly types
        for ipoly in range(self.object.npoly):
            if self.object.type[ipoly] > 0:
                if is_clockwise(self.object.poly[ipoly]):
                    self.object.type[ipoly] = 1
                else:
                    self.object.type[ipoly] = 2
        return self.object.head

    def put_content(self, file_name, head=None):
        # ~~> all object should have a put_content method
        self.object.put_content(file_name, head)

    def get_areas(self, select=None):
        if select == None:
            select = range(self.npoly)
        # TODO: take sperical coordinate into account
        areas = np.zeros(len(select), dtype=np.float)
        for i, ipoly in enumerate(select):
            # ~~> compute the area only for closed polygons
            if self.object.type[ipoly] in [1, 2]:
                areas[i] = get_area(self.object.poly[ipoly])
        return areas

    def get_lenghts(self, select=None):
        if select == None:
            select = range(self.npoly)
        return np.zeros(len(select), dtype=float)

    def sort_by_areas(self):
        srt = np.sort(np.array(zip(np.argsort(self.get_areas())[::-1],
                                   np.arange(self.npoly)),
                               dtype=[('s', int), ('i', int)]),
                      order='s')
        for i in range(len(srt)):
            if srt['s'][i] > srt['i'][i]:
                self.object.poly.insert(srt['i'][i],
                                        self.object.poly.pop(srt['s'][i]))
                self.object.type.insert(srt['i'][i],
                                        self.object.type.pop(srt['s'][i]))
                self.object.vals.insert(srt['i'][i],
                                        self.object.vals.pop(srt['s'][i]))

    def make_anti_clockwise(self, select=None):
        if select is None:
            select = range(self.npoly)
        for ipoly in select:
            # ~~> make anti-clockwise only closed clockwise polygons
            if self.object.type[ipoly] in [1, 2]:
                self.object.poly[ipoly], self.object.vals[ipoly] = \
                         make_anti_clockwise(self.object.poly[ipoly],
                                             self.object.vals[ipoly])
                self.object.type[ipoly] = 2

    def make_clockwise(self, select=None):
        if select is None:
            select = range(self.npoly)
        for ipoly in select:
            # ~~> make clockwise only closed anti-clockwise polygons
            if self.object.type[ipoly] in [1, 2]:
                self.object.poly[ipoly], self.object.vals[ipoly] = \
                        make_clockwise(self.object.poly[ipoly],
                                       self.object.vals[ipoly])
                self.object.type[ipoly] = 1

    def smooth_subdivise(self, select=None, weigth=0.5):
        if select is None:
            select = range(self.npoly)
        for ipoly in select:
            # ~~> make clockwise only closed anti-clockwise polygons
            self.object.poly[ipoly], self.object.vals[ipoly], _ = \
                      smooth_subdivise(self.object.poly[ipoly],
                                       self.object.vals[ipoly],
                                       self.object.type[ipoly],
                                       weigth)

    def subsample_distance(self, select=None, distance=1000.0):
        if select is None:
            select = range(self.npoly)
        for ipoly in select:
            # ~~> make clockwise only closed anti-clockwise polygons
            self.object.poly[ipoly], self.object.vals[ipoly], _ = \
                    subsample_distance(self.object.poly[ipoly],
                                       self.object.vals[ipoly],
                                       self.object.type[ipoly],
                                       distance)

    def subsample_angle(self, select=None, angle=15.0):
        if select is None:
            select = range(self.npoly)
        for ipoly in select:
            # ~~> make clockwise only closed anti-clockwise polygons
            self.object.poly[ipoly], self.object.vals[ipoly], _ = \
                      subsample_angle(self.object.poly[ipoly],
                                      self.object.vals[ipoly],
                                      self.object.type[ipoly],
                                      angle)

    def sph2ll(self, t1):
        (long0, lat0) = t1
        radius = 6371000.
        long0 = np.deg2rad(float(long0))
        lat0 = np.deg2rad(float(lat0))
        const = np.tan(lat0/2. + np.pi/4.)
        for poly in self.object.poly:
            for ipoly in range(len(poly)):
                poly[ipoly][0] = np.rad2deg(poly[ipoly][0]/radius + long0)
                poly[ipoly][1] = np.rad2deg\
                          (2.*np.arctan(const*np.exp(poly[ipoly][1]/radius)) \
                           - np.pi/2.)

    def ll2sph(self, t1):
        (long0, lat0) = t1
        radius = 6371000.
        long0 = np.deg2rad(float(long0))
        lat0 = np.deg2rad(float(lat0))
        const = np.tan(lat0/2. + np.pi/4.)
        for poly in self.object.poly:
            for ipoly in range(len(poly)):
                poly[ipoly][0] = radius * (np.deg2rad(poly[ipoly][0]) - long0)
                poly[ipoly][1] = radius * (np.log(\
                             np.tan(np.deg2rad(poly[ipoly][1])/2. + np.pi/4.)) \
                                                - np.log(const))

    def get_bbox(self):
        if len(self.object.poly) == 0:
            return 0., 0., 0., 0.
        xmin = xmax = self.object.poly[0][0][0]
        ymin = ymax = self.object.poly[0][0][1]
        for poly in self.object.poly:
            x_p, y_p = poly.T
            xmin = min(xmin, min(x_p))
            xmax = max(xmax, max(x_p))
            ymin = min(ymin, min(y_p))
            ymax = max(ymax, max(y_p))
        return xmin, ymin, xmax, ymax

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

__author__ = "Sebastien E. Bourban"
__date__ = "$15-Nov-2011 08:51:29$"

def main():
    """ Main function of polygons """

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Command line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nInterpreting command line options\n'+'~'*72+'\n')
    parser = ArgumentParser(\
        formatter_class=RawDescriptionHelpFormatter,
        description=('''\n
Testing ...
        '''))
    parser.add_argument("args", nargs='+')
    parser.add_argument("--sph2ll", dest="sph2ll", default=None,
                        help="convert from spherical to longitude-latitude")
    options = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads code name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ins_file = options.args[0]

    poly = Polygons()
    head = poly.parse_content(ins_file)

    poly.subsample_distance(distance=250.)
    poly.smooth_subdivise(weigth=0.7)
    poly.subsample_angle(angle=12.0)

    poly.put_content(options.args[1], head)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)

if __name__ == "__main__":
    main()
