#!/usr/bin/env python3
""" Handling shapefile to extract polygon """

from os import path

import data_manip.formats.shapefile as shapefile
from utils.exceptions import TelemacException


def read_shape_epsg(shape_file):
    """
    Extract a polygon from a shape file

    @param shape_file (string) Name of the shape file

    @return (str) the EPSG code

    """
    epsg = ''
    try:
        import fiona
    except ImportError:
        print('=> Information: fiona Python package is not available')
        print('   EPSG code cannot be recovered from Telemac')
        return epsg

    shp = fiona.open(shape_file, 'r')
    if 'init' in shp.crs:
        epsg = shp.crs['init']

    return epsg


def read_shape_data(shape_file, get_names=False):
    """
    Extract polygons or their names from a shape file

    @param shape_file (string) Name of the shape file
    @param get_names (boolean) If True returns the name of each shape

    @return (list) The list representing the polygon (list of 2-uple)
                   Or the list of the names of the polygon (if get_names=True)
    """
    isfiona = True
    try:
        import fiona
    except ImportError:
        print('=> Information: fiona Python package is not available')
        print('   Telemac shapefile implementation will be used instead')
        isfiona = False

    if isfiona:
        ext = shape_file.split('.')[1]
        shp = fiona.open(shape_file, 'r')
        if ext in ('shp', 'geojson', 'json'):
            if get_names:
                names = []
                for shp_obj in shp:
                    obj = shp_obj['properties']
                    if 'Name' in obj:
                        names.append(obj['Name'])
                    else:
                        names.append(None)
                return names

            res = []
            for shp_obj in shp:
                obj = shp_obj['geometry']
                if obj is None:
                    continue
                res.append(obj['coordinates'])
            return res
        if ext in ('gpkg', 'dxf'):
            if get_names:
                names = []
                for obj in shp:
                    if 'Name' in obj['properties']:
                        names.append(obj['properties']['Name'])
                    else:
                        names.append(None)
                return names

            res = []
            for obj in shp:
                if obj['geometry'] is None:
                    continue
                res.append(obj['geometry']['coordinates'])
            return res

    if get_names:
        root, _ = path.splitext(shape_file)
        shape = shapefile.Reader(shape_file, dbf=root+'.dbf')
    else:
        shape = shapefile.Reader(shape_file)

    if get_names:
        names = []
        for record in shape.iterRecords():
            if record != [None]:
                names.append(record[1])
        return names

    res = []
    # first feature of the shapefile
    for feature in shape.shapeRecords():
        first = feature.shape.__geo_interface__
        if first['type'] not in [("LineString"), ("Polygon"),
                                 ("Point")]:
            raise TelemacException(
                "No linestring or Polygon or Point in shapefile\n"
                "type: {} not handled".format(first['type']))

        if first['type'] == ("LineString"):
            res.append(list(first['coordinates']))
        elif first['type'] == ("Polygon"):
            res.append(list(first['coordinates'][0][:-1]))
        elif first['type'] == ("Point"):
            res.append(list(first['coordinates']))
    return res

def read_poly_txt(filename):
    """
    Read polygons stored in txt format (Telemac in_poly)

    @param filename(str) name of the file to read
    """
    with open(filename, 'r') as ffile:
        line = ffile.readline().strip()
        npolys = int(line)

        polys = []
        line = ffile.readline().strip()
        size_poly = [int(val) for val in line.split()]

        for i in range(npolys):
            poly = []
            for _ in range(size_poly[i]):
                line = ffile.readline().strip()
                point = [float(val) for val in line.split()]
                poly.append(point)
            polys.append(poly)

    return polys

def write_poly_txt(polys, filename):
    """
    Write some polygones in a txt file

    @param polys(list) list of the polygones to write
    @param filename(str) name of the file to save

    """
    with open(filename, 'w') as ffile:
        ffile.write("{}\n".format(len(polys)))
        for poly in polys:
            ffile.write("{} ".format(len(poly)))
        ffile.write("\n")
        for poly in polys:
            for point in poly:
                ffile.write("{} {} \n".format(point[0], point[1]))

#TODO: Have a look at data_manip.extraction.parser_kenue that seems to be
#      handling i2s and i3s files better
def write_poly_i2s(polys, filename, poly_val=None):
    """
    Write some polygones in a txt file

    @param polys (list) list of the polygones to write
    @param filename (str) name of the file to save
    @param poly_val (list) list of value for each polygon (if None will set to
    42.)
    """
    npoly = len(polys)
    if poly_val is None:
        poly_val = [42.]*npoly
    else:
        nval = len(poly_val)
        if nval != npoly:
            raise TelemacException(
                "Wrong lenght for poly_val is {} should be {}"
                .format(nval, npoly))

    with open(filename, 'w') as ffile:
        ffile.write(":EndHeader\n")
        for i, poly in enumerate(polys):
            ffile.write("{} {} \n".format(len(poly), poly_val[i]))
            for point in poly:
                ffile.write("{} {} \n".format(point[0], point[1]))
            ffile.write("{} {} \n".format(poly[0][0], poly[0][1]))

def write_shape(file_name, polys, poly_names=None, output_driver='ESRI Shapefile', crs=None):
    """
    Write list of polygons in Shape format

    @param file_name (str) Name of the shape file to write
    @param polys (list) list of list of points
    @param poly_names (list) List of name for each polygon
    @param output_driver (str) Format of output file
    @param crs (str/dict) Projection for the file
    """

    from collections import OrderedDict
    import fiona

    if poly_names is None:
        poly_names = ['polygon_{}'.format(i) for i in range(len(polys))]

    schema = {'geometry': 'Polygon',
              'properties': OrderedDict([('name', 'str')])
             }
    if crs is not None:
        shp = fiona.open('contour.shp', 'w', driver=output_driver,
                         schema=schema)
    else:
        shp = fiona.open('contour.shp', 'w', driver=output_driver,
                         schema=schema, crs=crs)

    for i, poly in enumerate(polys):
        shp_poly = {'geometry': {'type': 'Polygon',
                                 'coordinates': poly},
                    'properties': OrderedDict([('name', poly_names[i])])
                   }

        shp.write(shp_poly)

    shp.close()
