""" Function for extraction of contour of mesh """
import numpy as np
import matplotlib.path as mpltPath
from data_manip.extraction.telemac_file import TelemacFile

def detecting_boundaries(connectivity_table, npoints):
    """
    Detect boundary points with the connectivity table

    @param connectivity table (np.array of shape (npoints,3)) contains indexes
    nodes of each element of the mesh
    @param npoints (integer) number of nodes of the mesh
    points in the polyline
    @return connectivity_bnd_elt (np.array of shape (number of boundary
    element, 3)) a numpy array which contains indexes of nodes constituting
    boundary elements (left, right and center nodes)
    @return bnd_points (np.array of shape (number of boundary points, 1)) a
    numpy array wich contains indexes of all boundary points of the mesh
    """
    buf = []
    connectivity_bnd_elt = []

    div = (npoints - npoints % 10) // 10
    if div == 0:
        div = 1

    if np.__version__ < '1.19.0':
        import collections

        for i in range(npoints):
            if i % div == 0:
                print(i // div * 10, '%')

            connectivity_rows = collections.Counter(
                connectivity_table[
                    np.where((connectivity_table[:, 0] == i)
                             + (connectivity_table[:, 1] == i)
                             + (connectivity_table[:, 2] == i))
                    ].flatten())


            temp = np.array(
                list(connectivity_rows.keys()))[np.where(
                    np.array(list(connectivity_rows.values())) == 1)].tolist()
            buf.extend(temp)

            if temp != []:
                #to handle overconstrained element
                if i in temp:
                    temp.remove(i)
                temp.append(i)
                connectivity_bnd_elt.append(temp)

    else:
        for i in range(npoints):
            if i % div == 0:
                print(i // div * 10, '%')

            connectivity_rows = connectivity_table[
                np.where((connectivity_table[:, 0] == i)
                         + (connectivity_table[:, 1] == i)
                         + (connectivity_table[:, 2] == i))
                ].flatten()

            uniq, count = np.unique(connectivity_rows, return_counts=True)
            temp = uniq[np.where(count == 1)].tolist()
            buf.extend(temp)

            if temp != []:
                #to handle overconstrained element
                if i in temp:
                    temp.remove(i)
                temp.append(i)
                connectivity_bnd_elt.append(temp)


    buf = np.array(buf)
    connectivity_bnd_elt = np.array(connectivity_bnd_elt)

    bnd_points = np.unique(buf)

    return connectivity_bnd_elt, bnd_points

def detecting_boundaries_with_bnd_file(connectivity_bnd_table, bnd_points):
    """
    Identify boundaries when you have a boundary file
    """
    buf = []
    connectivity_bnd_elt = []


    div = (len(bnd_points) - len(bnd_points) % 10) // 10
    if div == 0:
        div = 1

    if np.__version__ < '1.19.0':
        import collections

        for i in bnd_points:
            if i % div == 0:
                print(i // div * 10, '%')

            connectivity_rows = collections.Counter(connectivity_bnd_table[
                np.where((connectivity_bnd_table[:, 0] == i)
                         + (connectivity_bnd_table[:, 1] == i)
                         + (connectivity_bnd_table[:, 2] == i))
                ].flatten())


            temp = np.array(
                list(connectivity_rows.keys()))[np.where(
                    np.array(list(connectivity_rows.values())) == 1)].tolist()

            if temp != []:
                #to handle overconstrained element
                if i in temp:
                    temp.remove(i)
                temp.append(i)
                connectivity_bnd_elt.append(temp)



    else:
        for j, i in enumerate(bnd_points):
            if j % div == 0:
                print(j // div * 10, '%')

            connectivity_rows = connectivity_bnd_table[
                np.where((connectivity_bnd_table[:, 0] == i)
                         + (connectivity_bnd_table[:, 1] == i)
                         + (connectivity_bnd_table[:, 2] == i))
                ].flatten()

            uniq, count = np.unique(connectivity_rows, return_counts=True)
            temp = uniq[np.where(count == 1)].tolist()
            buf.extend(temp)

            if temp != []:
                #to handle overconstrained element
                if i in temp:
                    temp.remove(i)
                temp.append(i)
                connectivity_bnd_elt.append(temp)


    connectivity_bnd_elt = np.array(connectivity_bnd_elt)

    return connectivity_bnd_elt

def get_first_point(tel, bnd_points):
    """
    Determine the southwest points among points given trough their indexes

    @param tel (TelemacFile) a Telemac file
    @param bnd_points (numpy array of shape (number of points)) contains indexes
    of points
    @return first_bnd_pt_index (integer) index of the southwest point
    """

    x_plus_y = tel.meshx[bnd_points] + tel.meshy[bnd_points]

    southwest_bnd_pts_index = \
        bnd_points[np.where(x_plus_y == x_plus_y.min())[0]]

    if southwest_bnd_pts_index.shape[0] == 1:
        first_bnd_pt_index = southwest_bnd_pts_index[0]
    else:
        first_bnd_pt_index = \
          southwest_bnd_pts_index[
                  np.where(tel.meshx[southwest_bnd_pts_index]
                      == tel.meshx[southwest_bnd_pts_index].min())[0][0]]

    return first_bnd_pt_index


def sorting_boundary(tel, first_pt_idx, connectivity_bnd_elt, clockwise=True):
    """
    Sort list of boundary nodes in clockwise or anticlockwise
    Connectivity table of the boundary elements is necessary

    @param tel (TelemacFile) a Telemac file
    @param first_bnd_pt_index (integer) index of the southwest point
    @param connectivity_bnd_elt (np.array of shape (number of boundary element,
    3)) a numpy array which contains indexes of nodes constituting boundary
    elements (left, right and center nodes)
    @param clockwise (boolean) direction of the sorting (optional, default is
    True)
    @return bnd (list) list of point indexes sorting clockwise or antclockwise
    """

    bnd = [first_pt_idx]
    neighbours_idx = connectivity_bnd_elt[
        np.where((connectivity_bnd_elt[:, 2] == bnd[-1]))][0][:-1]
    if clockwise:
        if tel.meshx[neighbours_idx[0]] > tel.meshx[neighbours_idx[1]]:
            bnd.append(neighbours_idx[0])
        else:
            bnd.append(neighbours_idx[1])

    else:
        if tel.meshx[neighbours_idx[0]] < tel.meshx[neighbours_idx[1]]:
            bnd.append(neighbours_idx[0])
        else:
            bnd.append(neighbours_idx[1])

    prev_pt = bnd[-2]
    while bnd[-1] != first_pt_idx:
        neighbours_idx = connectivity_bnd_elt[
            np.where((connectivity_bnd_elt[:, 2] == bnd[-1]))][0][:-1]
        bnd.append(np.delete(neighbours_idx,
                             np.where((neighbours_idx == prev_pt)
                                      | (neighbours_idx == prev_pt)))[0])

        prev_pt = bnd[-2]

    return bnd

def sorting_boundaries(tel, bnd_points, connectivity_bnd_elt):
    """
    Sort boundaries in the case there is islands and detect
    the case of separated domains

    @param tel (TelemacFile) a Telemac file
    @param bnd_points (np.array of shape (number of boundary points, 1)) a numpy
    array wich contains indexes of all boundary points of the mesh
    @return connectivity_bnd_elt (numpy array of shape (number of boundary
    element, 3)) a numpy array which contains indexes of nodes constituting
    boundary elements (left, right and center nodes)
    @return boundaries (list) list of boundaries, external and islands if
    existing, which are list of point indexes (clockwise for external boundary,
    anticlockwise for islands)
    @return left_pts_idx_outside (numpy array) point indexes of boundary points
    not inside current external boundary analyzed
    @return left_bnd_elt_outside (numpy array of shape (number of boundary
    element, 3)) a numpy array which contains indexes of nodes constituting
    boundary elements (left, right and center nodes) not inside current
    external boundary analyzed
    """

    boundaries = []

    first_pt_idx = get_first_point(tel, bnd_points)
    boundaries.append(sorting_boundary(tel, first_pt_idx, connectivity_bnd_elt))
    if len(boundaries[0]) - 1 == len(bnd_points):
        print("No islands")
        return boundaries, np.array([]), np.array([])

    left_bnd_elt = connectivity_bnd_elt
    j = 0

    poly = np.column_stack((tel.meshx[boundaries[0]], tel.meshy[boundaries[0]]))
    path = mpltPath.Path(poly)
    for i in range(len(boundaries[j])):
        left_bnd_elt = np.delete(left_bnd_elt, np.where(
            (left_bnd_elt[:, 0] == boundaries[0][i]) |
            (left_bnd_elt[:, 1] == boundaries[0][i]) |
            (left_bnd_elt[:, 2] == boundaries[0][i])), 0)


    left_pts_idx = left_bnd_elt[:, 2]
    left_pts = np.column_stack((tel.meshx[left_pts_idx],
                                tel.meshy[left_pts_idx]))

    left_pts_idx_inside = left_pts_idx[path.contains_points(left_pts)]

    left_pts_idx_outside = left_pts_idx[
        np.invert(path.contains_points(left_pts))]
    if left_pts_idx_outside.size > 0:
        left_bnd_elt_outside = left_bnd_elt[
            np.where(left_bnd_elt[:, 2] == left_pts_idx_outside)]
    else:
        left_bnd_elt_outside = np.array([])

    if left_pts_idx_inside.size > 0:
        left_bnd_elt_inside = left_bnd_elt[
            np.where(left_bnd_elt[:, 2] == left_pts_idx_inside)]
        print("Islands")
        while left_bnd_elt_inside.shape[0] != 0:
            j += 1

            #islands case
            next_first = np.where(left_pts_idx_inside ==
                                  get_first_point(tel,
                                                  left_bnd_elt_inside[:, 2]))
            boundaries.append(sorting_boundary(
                tel,
                left_pts_idx_inside[next_first][0],
                left_bnd_elt_inside,
                clockwise=False))


            for i in range(len(boundaries[j])):
                left_bnd_elt_inside = np.delete(left_bnd_elt_inside, np.where(
                    (left_bnd_elt_inside[:, 0] == boundaries[j][i]) |
                    (left_bnd_elt_inside[:, 1] == boundaries[j][i]) |
                    (left_bnd_elt_inside[:, 2] == boundaries[j][i])), 0)

            left_pts_idx_inside = left_bnd_elt_inside[:, 2]




        return boundaries, left_pts_idx_outside, left_bnd_elt_outside

    print("No Islands")

    return boundaries, left_pts_idx_outside, left_bnd_elt_outside


def write_gis_file(polygons, \
                   file_name='contour.shp', \
                   output_driver='ESRI Shapefile', \
                   crs=None):
    """
    Write polygon GIS file (multiple polygon is allowed)

    @param polygon (list of lists of tuples) list of polygons which are
    described by their vertex coordinates (external vertices are clockwise
    sorting and holes are anticlockwise sorting)
    @param file_name (string) name of the written file (can contain a path with
    os.sep)
    @param crs_given (boolean) optional parameter to indicate the presence of a
    coordinate reference system (default: False)
    @param output_driver (string) define the format of the GIS file used by
    fiona.  To know available drivers, use fiona.available_drivers attribute of
    Fiona library
    @param crs (string) EPSG (European Petroleum Survey Group) code of the
    coordinate reference system
    """
    from collections import OrderedDict
    import fiona
    schema = {'geometry': 'Polygon',
              'properties': OrderedDict([('name', 'str')])
             }

    dict_polygons = []

    for i, polygon in enumerate(polygons):

        dict_poly = {
            'geometry': {
                'type': 'Polygon',
                'coordinates': polygon},
            'properties': OrderedDict([('name', 'polygon_' + str(i))])
            }

        dict_polygons.append(dict_poly)

    with fiona.open(file_name, 'w', driver=output_driver, crs=crs,
                    schema=schema) as shp:
        for dict_poly in dict_polygons:
            shp.write(dict_poly)

def extract_contour(tel_file_path, \
                    bnd_file=None, ):
    """
    Generic function for extraction of contour from a mesh (with our without
    boundary file)

    @param tel_file_path (str) Path to the mesh file
    @param bnd_file (str) Path to the boundary file

    @returns (list) List of polygons
    """

    domains = []

    if bnd_file is None:
        tel = TelemacFile(tel_file_path)
        print("Number of nodes: ", tel.npoin2)
        connectivity_bnd_elt, bnd_points =\
            detecting_boundaries(tel.ikle2, tel.npoin2)
    else:
        tel = TelemacFile(tel_file_path, bnd_file=bnd_file)
        print("Number of nodes: ", tel.npoin2)
        bnd_points = tel.get_bnd_numbering()
        connectivity_bnd_elt = detecting_boundaries_with_bnd_file(tel.ikle2,
                                                                  bnd_points)

        if len(bnd_points) == 0:
            print("Boundary file seems to not correspond to the mesh \n")
            print("Process like if there was no boundary file")
            connectivity_bnd_elt, bnd_points =\
                    detecting_boundaries(tel.ikle2, tel.npoin2)

    print("Number of boundary nodes: ", len(bnd_points))


    print("Domain 1")
    boundaries_ordered, bnd_idx_left, bnd_elt_left = \
        sorting_boundaries(tel, bnd_points, connectivity_bnd_elt)
    domains.append(boundaries_ordered)

    i = 1
    while bnd_elt_left.size != 0:
        i += 1
        print("Domain " + str(i))

        boundaries_ordered, bnd_idx_left, bnd_elt_left = \
            sorting_boundaries(tel, bnd_idx_left, bnd_elt_left)
        domains.append(boundaries_ordered)

    domains_bnd = []
    for domain in domains:
        poly_bnd = []
        for bnd in domain:
            coord_bnd = []
            for x, y in zip(tel.meshx[bnd], tel.meshy[bnd]):
                coord_bnd.append((x, y))
            poly_bnd.append(coord_bnd)
        domains_bnd.append(poly_bnd)

    return domains_bnd

def main():
    from utils.files import recursive_glob
    import os

    modules = ["telemac2d", "telemac3d", "sisyphe", "gaia", "tomawac"]

    for module in modules:
        med_files = recursive_glob(
            os.path.join(os.environ['HOMETEL'], 'examples', module),
            '*.med')
        slf_files = recursive_glob(
            os.path.join(os.environ['HOMETEL'], 'examples', module),
            '*.slf')

        geo_files = med_files + slf_files

        folder = os.path.join(".", module)

        if not os.path.isdir(folder):
            os.mkdir(folder)

        exception_prefix = ["f2d", "r2d", "t2d", "ini", "pre", "fsp",
                            "spe", "f3d", "r3d", "t3d"]

        for geo in geo_files:
            prefix = geo.split(os.sep)[-1].split('.')[-2]
            mesh_format = geo.split(os.sep)[-1].split('.')[-1]
            geo_name = prefix + '_' + mesh_format + '.shp'
            if prefix == "fom_spe_friction" or \
                prefix[0:3] in exception_prefix:
                continue

            geo_folder = os.sep.join(geo.split(os.sep)[:-1])
            geo_name = os.path.join(folder, geo_name)
            cli_file = os.path.join(geo_folder, prefix + ".cli")
            if os.path.isfile(cli_file):
                bnd_file = cli_file
            else:
                bnd_file = None
            if not os.path.isfile(geo_name):
                print(geo, bnd_file)
                domains_bnd = extract_contour(geo, \
                                bnd_file=bnd_file)

                write_gis_file(domains_bnd, file_name=geo_name)

                print("\n")

if __name__ == '__main__':
    main()
