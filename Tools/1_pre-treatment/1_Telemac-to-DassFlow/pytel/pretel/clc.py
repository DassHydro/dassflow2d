"""
Contains script to load CorineLandCover information
"""
from os import path, environ
import re
import numpy as np
from utils.exceptions import TelemacException
from utils.polygon import points_in_poly
from data_manip.extraction.telemac_file import TelemacFile

def add_clc_in_file(geo_file, varname, year, res_file=None):
    """
    Function to extract clc and write it into a file

    @param geo_file (str) Name of the input file
    @param res_file (str) Name of the output file
    """
    if not path.exists(geo_file):
        raise TelemacException("File not found: {}".format(geo_file))
    # Output file
    if res_file is not None and res_file != geo_file:
        if path.exists(res_file):
            raise TelemacException(
                "Output file already exists: {}".format(res_file))
        res = TelemacFile(res_file, access='w')
        geo = TelemacFile(geo_file, access='r')
    else:
        res = TelemacFile(geo_file, access='rw')
        geo = res

    # Extracting clc data
    # Xml file
    xml_file = 'clc_{}.xml'.format(year)
    #
    print(" "*6+" ~> Extraction wfs data into xml file")
    x_min = geo.meshx.min()
    y_min = geo.meshy.min()
    x_max = geo.meshx.max()
    y_max = geo.meshy.max()
    wfs2xml(xml_file, year, (x_min, y_min, x_max, y_max))

    print(" "*6+" ~> Building clc structure")
    struct = build_clc_struct(xml_file, year)

    print(" "*6+" ~> Building correspondace table")
    corresp = read_corresp_table(year)

    res.read(geo)

    print(" "*6+" ~> Writting data in mesh")
    # Computing friction coefficient on each mesh point
    friction = set_friction_clc_on_mesh(res.meshx, res.meshy, struct, corresp)

    ## Adding new variable FRICTION
    res.add_variable(varname, '')

    # Setting value to friction avlue
    res.add_data_value(varname, 0, friction)

    # writting new file
    res.write()

    res.close()


def wfs2xml(xml_file, year, bbox):
    """
    Load from a wfs flux the information and write them in a xml file

    @param xml_file (str) File to write
    @param bbox (4-upple) Coordinates of (bottom left point x,  and y, upper
    right point,  x and y)
    """
    typename = get_typename(year)

    from owslib.wfs import WebFeatureService

    wfs11 = WebFeatureService('http://wxs.ign.fr/corinelandcover/geoportail/' +
                              'wfs?service=WFS&request=GetCapabilities')
    response = wfs11.getfeature(typename=[typename], bbox=bbox)

    out = open(xml_file, 'wb')
    data = response.read()
    if isinstance(data, str):
        data = bytes(data, 'UTF-8')
    out.write(data)
    out.close()


def get_typename(year):
    """
    Return typename in the wfs from the year

    @param year (int) Year to extract
    """
    if year == 2012:
        typename = "LANDCOVER.CLC12_FR:clc12_fr"
    else:
        raise Exception("Unknown year: {}\nYear handled: 2012".format(year))

    return typename


def build_clc_struct(xml_file, year):
    """
    Read the xml file and extract information

    @param xml_file (str) File to write
    @param year (int) Year to extract
    """
    import xml.etree.ElementTree as et
    tree = et.parse(xml_file)
    root = tree.getroot()

    tag, _ = get_typename(year).split(':')

    struct = {}

    for child in root:
        if child.tag == '{http://www.opengis.net/gml}featureMember':
            zone = child[0]
            for zone_elem in zone:
                # Here we are in a new strickler zone
                # print(zone_elem.tag)
                if zone_elem.tag == \
                        '{http://wxs.ign.fr/datastore/'+tag+'}code_12':
                    # Saving code for that zone
                    code = zone_elem.text
                    if code not in struct:
                        struct[code] = []
                    # add a key to store apart the inner
                    # boundaries of a code (MS)
                    if code+'inner' not in struct:
                        # print(code+'inner')
                        struct[code+'inner'] = []
                    # print("code: ", code)
                zone_tag = '{http://wxs.ign.fr/datastore/'+tag+'}the_geom'
                if zone_elem.tag == zone_tag:
                    # Adding the polygons we have in there
                    polys_item = zone_elem[0][0][0]
                    for poly_item in polys_item:
                        poly_txt = poly_item[0][0].text
                        poly = np.array([[float(i) for i in elem.split(',')]
                                         for elem in poly_txt.split(' ')])
                        # Store outer boundaries and inner boundaries in
                        # different keys (MS)
                        if poly_item.tag == \
                            '{http://www.opengis.net/gml}outerBoundaryIs':
                            struct[code].append(poly)
                        else:
                            struct[code+'inner'].append(poly)
    return struct

def read_corresp_table(year, table_name=None):
    """
    Read the corresponding file

    @param year (int) Year to extract
    @parma table (str) Overwrites the table used
    """
    if year == 2012:
        tab_file = path.join(environ['HOMETEL'], 'notebooks', 'data',
                             'def_strickler_table_12.txt')
    else:
        raise Exception("Year not handled")

    if table_name is not None:
        tab_file = table_name

    corresp = {}
    with open(tab_file, 'r') as f:
        f.readline()
        for line in f.readlines():
            match = re.match(r'"(?P<name>[^"]+)"\s+(?P<value>[^ ]+)\s+(?P<ref>[^ ]+)\s+(?P<code>[^ ]+)', line)
            corresp[int(match.group('code'))] = {\
                    'name':match.group('name'),
                    'ref':match.group('ref'),
                    'value':match.group('value')}

    return corresp


def plot_clc(clc_struct, corresp, fig_name=None, xmin=None, xmax=None, \
            ymin=None, ymax=None):
    """
    Plotting all the polyline from the clc structure

    @param clc_struct (dico) Structure of the clc data
    """
    import matplotlib.pyplot as plt

    _, ax = plt.subplots()

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']*7

    for j, code in enumerate(clc_struct):
        if 'inner' in code:
            continue
        color = colors[j]
        for i, poly in enumerate(clc_struct[code]):
            value = corresp[int(code)]['value']
            if i == 0:
                ax.plot(poly[:, 0], poly[:, 1], color=color,
                        label='code {} ({})'.format(code, value))
            else:
                ax.plot(poly[:, 0], poly[:, 1], color=color)

    chart_box = ax.get_position()
    ax.set_position([chart_box.x0, chart_box.y0,
                     chart_box.width*0.6, chart_box.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.45, 1.0),
              shadow=True, ncol=1)

    if xmin is not None:
        plt.xlim(xmin, xmax)
    if ymin is not None:
        plt.ylim(ymin, ymax)
    if fig_name is not None:
        print(" ~> Plotting {}".format(fig_name))
        plt.savefig(fig_name)
    else:
        plt.show()

    plt.close('all')

def set_friction_clc_on_mesh(meshx, meshy, clc_struct, corresp):
    """
    Compute friction for each point of mesh_file according to clc into file

    @param meshx (np.array) X coordinates pf mesh points
    @param meshy (np.array) Y coordinates pf mesh points
    @param clc_strucut (dict) Structure for the clc
    @param corresp (dict) corresponding table between clc code and value
    """

    npoin = meshx.shape[0]
    friction = np.zeros((npoin), dtype=np.float64)
    points = np.zeros((npoin, 2), dtype=np.float64)
    points[:, 0] = meshx
    points[:, 1] = meshy
    for code in clc_struct:
        # Test to treat only "code" keys and not "code"+inner keys (MS)
        if code[-1] != 'r':

            for poly in clc_struct[code]:

                points_in = points_in_poly(points, poly)
                friction[points_in] = corresp[int(code)]['value']

    return friction
