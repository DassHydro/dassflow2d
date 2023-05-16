"""
tel2tom and tom2tel
"""
from os import path, remove
import numpy as np
import matplotlib.tri as mtri
from data_manip.extraction.telemac_file import TelemacFile
from utils.polygon import points_in_poly

def interp_triangle_area(xtri, ytri, x, y):
    """
     performs barycentrix interpolation (similar to bilinear)
    xTri: x-coordinates of triangle (Nx3)
    yTri: y-coordinates of triangle (Nx3)
     x = x coordinates for interpolation (Nx1)
     y = y coordinates for interpolation (Nx1)
    """
    a_1 = 0.5*np.abs((x-xtri[:, 2])*(ytri[:, 1]-y) -
                     (x-xtri[:, 1]) * (ytri[:, 2]-y))
    a_2 = 0.5*np.abs((xtri[:, 0]-xtri[:, 2])*(y-ytri[:, 0]) -
                     (xtri[:, 0]-x) * (ytri[:, 2]-ytri[:, 0]))
    a_3 = 0.5*np.abs((xtri[:, 0]-x)*(ytri[:, 1]-ytri[:, 0]) -
                     (xtri[:, 0]-xtri[:, 1]) * (y-ytri[:, 0]))
    a_a = a_1 + a_2 + a_3

    return np.vstack((a_1/a_a, a_2/a_a, a_3/a_a)).T


def interp_triangle_prepare(connection, xtri, ytri, x, y,
                            tree=None, extrap=False):
    """
     This function determines all preprocessing for trainagle interpolation

     sctInterp = interpTrianglePrepare(connection,xTri,yTri,x,y)

    INPUT
       - xTri: x-coordinates of triangle (Mx1)
       - yTri: y-coordinates of triangle (Mx1)
       - connection: coordinate number for interpolation (Kx3)
       - x = x coordinates for interpolation (Nx1)
       - y = y coordinates for interpolation (Nx1)
       - haveWaitBar: logical to determine wethre a waitbar is shown(optional)
       - extrap: logical to determine whether nearest neighbour interpolation
                 is done (optional)

    OUTPUT
     CoordIndex of the points on the triangle (Nx3) needed for interpolation
     interpCoef (coefficient) used to multiply
     mask: logical value of points
    """
    # determine in which triangle the points are

    my_tri = mtri.Triangulation(xtri, ytri, connection)
    tri_finder = my_tri.get_trifinder()

    intri = tri_finder(x, y)

    # Initialising
    coordx = np.empty((len(x), 3), dtype=np.float64)
    coordy = np.empty((len(x), 3), dtype=np.float64)
    coord_index = np.empty((len(x), 3), dtype=np.float64)

    coordx[:] = np.nan
    coordy[:] = np.nan
    coord_index[:] = np.nan

    # Excluding out of mesh points
    mask = intri != -1
    nomask = np.logical_not(mask)

    for i in range(3):
        coord_index[mask, i] = connection[intri[mask], i]
        coordx[mask, i] = xtri[connection[intri[mask], i]]
        coordy[mask, i] = ytri[connection[intri[mask], i]]

    interp_coef = interp_triangle_area(coordx, coordy, x, y)

    if extrap:
        extrap_x = np.array([[xx, yy] for xx, yy in zip(x, y)])
        _, index = tree.query(extrap_x[nomask, :])
        coord_index[nomask, :] = [0, -1, -1]
        coord_index[nomask, 0] = index
        interp_coef[nomask, :] = [1., 0., 0.]

    return coord_index, interp_coef


def connect_tel2tom(tel_file, tom_file, contour_tom=None, contour_tel=None):
    """
    connectTelTom(telFile,tomFile,contourTom,contourTel)

    INPUT:
    - telFile, tomFile: selafin files with the meshes of telemac
    and tomawac
    - contourTel, contourTom: i2s files with contours that should
    not be taken into account
    """
    t2d = TelemacFile(tel_file)
    tom = TelemacFile(tom_file)

    x_t2d = t2d.meshx
    y_t2d = t2d.meshy
    x_tom = tom.meshx
    y_tom = tom.meshy
    xy_t2d = np.vstack((x_t2d, y_t2d)).T
    xy_tom = np.vstack((x_tom, y_tom)).T

    # tel2tom
    print("  ~> Building tel2tom")
    if contour_tel is not None:
        # Select points within contour
        mask = points_in_poly(xy_t2d, contour_tel)
    else:
        mask = np.zeros(len(x_tom), dtype=np.bool)
        mask[:] = True

    res_file = path.splitext(tom_file)[0]+'-tel2tom.slf'
    root, ext = path.splitext(tom_file)
    res_file = root+'-tel2tom'+ext
    if path.exists(res_file):
        remove(res_file)
    res = TelemacFile(res_file, access='w')
    connect_send_recv('TEL2TOM', mask, t2d, tom, res)
    res.close()
    print("  - Created {}".format(res_file))

    # tom2tel
    print("  ~> Building tom2tel")
    if contour_tom is not None:
        # Select points within contour
        mask = points_in_poly(xy_tom, contour_tom)
    else:
        mask = np.zeros(len(x_t2d), dtype=np.bool)
        mask[:] = True

    root, ext = path.splitext(tel_file)
    res_file = root+'-tom2tel'+ext
    if path.exists(res_file):
        remove(res_file)
    res = TelemacFile(res_file, access='w')
    connect_send_recv('TOM2TEL', mask, tom, t2d, res)
    res.close()
    print("  - Created {}".format(res_file))


    t2d.close()
    tom.close()


def connect_send_recv(var_name, mask, send, recv, out):
    """
     computed connectivity list between telemac and tomowac

    connectSendRecv(varName,mask,sctSend,sctRecv,outFile)

     INPUT:
     - varName: name of variables to create. can be either TEL2TOm
     or TOM2TEL
     - mask: list of points not to include
     - sctSend: - strcuture of mesh for sender using telheadr
     - sctRecv: - strcuture of mesh for receiver using telheadr
     - outFile: name of file to make

     IT COMPUTES THE CONNECTIVITY LIST BETWEEN
     TELEMAC2D AND TOMAWAC SELAFIN FILES

     --------------------------------------------------------------
     The new TOMAWAC selafin file TOM_interp.slf created, includes:
     --------------------------------------------------------------
     ->TEL2TOM:
     CLOSEST INDEX OF THE TELEMAC2D MESH ONTO
     THE TOMAWAC GRID NODES
     ->TEL2TOM01,TEL2TOM02,TEL2TOM03:
     NEAREST NEIGHBOR INDEX (THE OTHER TWO VARIABLES ARE ZERO)
     ->TEL2TOMWTS01,TEL2TOMWTS02,TEL2TOMWTS03:
     LINEAR INTERPOLATION COEFFICIENTS
    """
    send.set_kd_tree()
    coord_index, interp_coef = interp_triangle_prepare(\
            send.ikle2, send.meshx, send.meshy,
            recv.meshx[mask], recv.meshy[mask],
            tree=send.tree, extrap=True)

    out.read(recv)

    for i in range(1, 4):
        varname = "{}{:02d}".format(var_name, i)
        out.add_variable(varname, '')
        out.add_data_value(varname, 0, coord_index[:, i-1] + 1)

        varname = "{}WTS{:02d}".format(var_name, i)
        out.add_variable(varname, '')
        out.add_data_value(varname, 0, interp_coef[:, i-1])

    out.write()
