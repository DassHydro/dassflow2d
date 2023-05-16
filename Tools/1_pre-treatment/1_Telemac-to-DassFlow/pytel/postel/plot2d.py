"""
Contains functions to plot scalar and vectors in 2d
"""
from collections import OrderedDict
from data_manip.formats.regular_grid import interpolate_on_grid
from postel.deco_cbar import deco_cbar
import numpy as np
import matplotlib.cm as cm
import matplotlib.tri as mtri
import matplotlib as mpl
import matplotlib.pyplot as plt


def set_frame_limits(mesh, margin=0.02):
    """
    Computes frame limits based on mesh x,y coordonates

    @param mesh (matplotlib.tri.Triangulation or tuple)
    triangular mesh or regular x,y grid
    @param margin (float) percentage used to set frame margin (default: 0.02)
    """
    if isinstance(mesh, mtri.triangulation.Triangulation):
        x, y = mesh.x, mesh.y
    elif isinstance(mesh, tuple):
        x, y = mesh
    else:
        raise ValueError("Unknown mesh type")

    xlen = np.max(x) - np.min(x)
    ylen = np.max(y) - np.min(y)
    xlim = [np.min(x) - margin*xlen, np.max(x) + margin*xlen]
    ylim = [np.min(y) - margin*ylen, np.max(y) + margin*ylen]

    return xlim, ylim


def set_extrema(data, vmin, vmax, eps=0.):
    """
    Set default min and max values for display range

    @param data (np.array) Data values
    @param vmin (float) Minimal value of data given
    @param vmax (float) Maximal value of data given
    @param eps (float) epsilon used to adjust extrema (default: 0.)
    """
    if vmin is None:
        vmin = np.min(data) - eps
    if vmax is None:
        vmax = np.max(data) + eps
    return vmin, vmax


def mask_triangles(tri, data, relation='leq', threshold=0.001):
    """
    Defines a boolean mask to mask triangles
    when data on triangles match criteria (relation to a threshold)

    @param tri (matplotlib.tri.Triangulation) triangular mesh
    @param data (numpy.array) scalar data used for masking
    @param threshold (float) threshold used for masking
    @param relation (str) how the masking is done relatively to the threshold,
    can be 'leq' (=<), 'geq' (=>), 'neq' (!=) or 'eq' (==), (default 'leq')
    """
    out = []
    for points in tri.triangles:
        i, j, k = points
        data_tri = (data[i] + data[j] + data[k])/3.
        if relation == 'geq':
            if data_tri >= threshold:
                out.append(True)
            else:
                out.append(False)
        elif relation == 'leq':
            if data_tri <= threshold:
                out.append(True)
            else:
                out.append(False)
        elif relation == 'eq':
            if data_tri == threshold:
                out.append(True)
            else:
                out.append(False)
        elif relation == 'neq':
            if data_tri != threshold:
                out.append(True)
            else:
                out.append(False)
        else:
            raise ValueError("Unknown relation, try 'leq', 'geq' or other")
    return out


def rotate(x, y, angle=0., center=[0., 0.]):
    """
    Clockwise rotation of a field (x,y)

    @param x (list) list of x-coordonate of the vectors
    @param y (list) list of y-coordonate of the vectors
    @param angle (float) angle of clockwise rotation in degree
    @param center (list) (default=[0., 0.]) center of rotation
    """
    rot00 = np.cos(np.pi*angle/180.)
    rot01 = np.sin(np.pi*angle/180.)
    xnew = rot00*(x - center[0]) + rot01*(y - center[1])
    ynew = -rot01*(x - center[0]) + rot00*(y - center[1])
    return xnew + center[0], ynew + center[1]


def plot2d_annotate_bnd(axe, tri, bnd_info,
                        colors=None, **kwargs):
    """
    Annotate boundary conditions

    @param axe (matplotlib.axes) matplotlib axes on which to draw
    @param tri (matplotlib.tri.Triangulation) triangular mesh
    @param colors (list) colors associated to boundary conditions
    (default plt.rcParams['axes.prop_cycle'].by_key()['color'])
    @param bnd_info (tuple) boundary conditions information (default None)
    @param kwargs (dict) rest of optional arguments given to plot
    """
    x, y = tri.x, tri.y
    nbor, lihbor, liubor, livbor, _ = bnd_info
    nbnd_poin = len(nbor)

    if colors is None:
        colors = mpl.pyplot.rcParams['axes.prop_cycle'].by_key()['color']

    bnd_types_dict = {'Closed boundaries/walls (2,2,2)': [2, 2, 2],
                      'Prescribed H (5,4,4)':            [5, 4, 4],
                      'Prescribed Q (4,5,5)':            [4, 5, 5],
                      'Prescribed Q and H (5,5,5)':      [5, 5, 5],
                      'Prescribed UV (4,6,6)':           [4, 6, 6],
                      'Prescribed UV and H (5,6,6)':     [5, 6, 6],
                      'Incident waves (1,1,1)':          [1, 1, 1],
                      'Custom (0,0,0)':                  [0, 0, 0],
                      'Free boundaries (4,4,4)':         [4, 4, 4]}
    label_bool = [False for i in range(len(bnd_types_dict))]

    labels = []
    handles = []
    for i in range(nbnd_poin):
        idx = 0
        bc_type = [lihbor[i], liubor[i], livbor[i]]
        colori = 'k'
        found = False
        for key in sorted(bnd_types_dict.keys()):
            item = bnd_types_dict[key]
            if bc_type == item:
                colori = colors[idx]
                # Adding to legends if never added
                if not label_bool[idx]:
                    label_bool[idx] = True
                    labels.append(key)
                    handles.append(mpl.lines.Line2D([0], [0],
                                   lw=0., marker='o', color=colori))
                found = True
                break
            idx += 1
        if not found:
            print(" ~> Could not find name for boundary type {}"
                  .format(bc_type))

        xi, yi = x[nbor[i]], y[nbor[i]]
        axe.plot(xi, yi, color=colori, **kwargs)

    axe.legend(handles, labels)


def plot2d_annotate_liq_bnd(axe, tri, liq_bnd_info, colors=None, **kwargs):
    """
    Annotate liquid boundaries

    @param axe (matplotlib.axes) matplotlib axes on which to draw
    @param tri (matplotlib.tri.Triangulation) triangular mesh
    @param colors (list) colors associated to boundary conditions
    (default plt.rcParams['axes.prop_cycle'].by_key()['color'])
    @param liq_bnd_info (tuple) boundary conditions information (default None)
    @param kwargs (dict) rest of optional arguments given to plot
    """
    x, y = tri.x, tri.y
    nbor, numliq = liq_bnd_info
    nbnd_poin = len(nbor)

    if colors is None:
        colors = mpl.pyplot.rcParams['axes.prop_cycle'].by_key()['color']

    nliq_bnd = np.max(numliq)

    # If we have at least one solid boundary
    if np.min(numliq) == 0:
        # Initialise labels with the solid boundary
        labels = ['Solid Boundary']
        handles = [mpl.lines.Line2D([0], [0], lw=0.,
                   marker='o', color=colors[0])]
    else:
        labels = []
        handles = []

    # Creating labels/handles for each liquid boundary (They are number from 1)
    for i in range(1, nliq_bnd+1):
        labels.append('Liquid boundary {}'.format(i))
        handles.append(mpl.lines.Line2D([0], [0],
                       lw=0., marker='o', color=colors[i]))

    # Plotting each boundary mapping it with its
    for i in range(nbnd_poin):
        x_i, y_i = x[nbor[i]], y[nbor[i]]
        axe.plot(x_i, y_i, color=colors[numliq[i]], **kwargs)

    axe.legend(handles, labels)


def plot2d_triangle_mesh(axe, tri, x_label='', y_label='',
                         color='k', linewidth=0.1,
                         **kwargs):
    """
    Plot a 2d representation of a triangle mesh

    @param axe (matplotlib.axes) matplotlib axes on which to draw
    @param tri (matplotlib.tri.Triangulation) triangular mesh
    @param x_label (string) Name of the x_label (default '')
    @param y_label (string) Name of the y_label (default '')
    @param color (str) mesh color (default 'k')
    @param linewidth (float) thickness of the plots line (default 0.1)
    @param kwargs (dict) rest of optional arguments given to triplot
    """
    axe.triplot(tri, color=color, linewidth=linewidth, **kwargs)
    xlim, ylim = set_frame_limits(tri)
    axe.set_xlim(xlim[0], xlim[1])
    axe.set_ylim(ylim[0], ylim[1])
    axe.set_xlabel(x_label)
    axe.set_ylabel(y_label)


def plot2d_image(axe, x_label='', y_label='',
                 image=None, image_file=None, **kwargs):
    """
    Plot a 2d image

    @param axe (matplotlib.axes) matplotlib axes on which to draw
    @param x_label (string) Name of the x_label (default '')
    @param y_label (string) Name of the y_label (default '')
    @param image (array) image to plot (default: None)
    @param image_file (str) image file to plot (default: None)
    @param kwargs (dict) rest of optional arguments given to imshow
    """
    assert image is not None or image_file is not None

    if image_file is not None:
        image = mpl.pyplot.imread(image_file)

    axe.imshow(image, **kwargs)
    axe.set_xlabel(x_label)
    axe.set_ylabel(y_label)


def plot2d_scalar_map(fig, axe, mesh, data,
                      x_label='', y_label='', data_name='data',
                      vmin=None, vmax=None, nv=None,
                      colorbar=True, cbar_properties=None, cbar_ticks=None,
                      cbar_ax=None, cbar_cax=None,
                      shading='gouraud', cmap_name='jet', **kwargs):
    """
    Plot the 2d map of scalar data on triangulation
    (matplotlib tripcolor)

    @param fig (matplotlib.figure) matplotlib figure structure
    @param axe (matplotlib.axes) matplotlib axes on which to draw
    @param mesh (matplotlib.tri.Triangulation or tuple)
    triangular mesh or regular x,y grid
    @param data (numpy.array) scalar data to plot
    @param x_label (string) Name of the x_label (default '')
    @param y_label (string) Name of the y_label (default '')
    @param data_name (str) name of the data to display
    @param vmin (float) Minimal value of data to plot
    @param vmax (float) Maximal value of data to plot
    @param nv (integer) Number of sample for colorbar range
    @param colorbar (bool) show colorbar (default: True)
    @param cbar_ticks (list) list of values where to show color bar ticks
    @param cbar_ax (Axes) Parent axes from which space for a new colorbar
    axes will be stolen. If a list of axes is given they will all be resized
    to make room for the colorbar axes.
    @param cbar_cax (Axes) Axes into which the colorbar will be drawn.
    @param cbar_properties (dict) list additional properties of the colorbar
    @param shading (str) shading used in tripcolor (default 'gouraud')
    @param cmap_name (str) colormap name (default 'jet')
    @param kwargs (dict) rest of optional arguments given to tripcolor
    """
    # Getting cmap object from name
    cmap = cm.get_cmap(name=cmap_name, lut=None)

    # Plotting scalar color map
    if isinstance(mesh, mtri.triangulation.Triangulation):
        img = axe.tripcolor(mesh, data, shading=shading,
                            cmap=cmap, vmin=vmin, vmax=vmax,
                            **kwargs)
    elif isinstance(mesh, tuple):
        x, y = mesh
        img = axe.pcolormesh(x, y, data, shading=shading,
                             cmap=cmap, vmin=vmin, vmax=vmax,
                             **kwargs)
    else:
        raise ValueError("Unknown mesh type")

    # Adding color bar if requested
    if colorbar:
        if nv is not None:
            vmin, vmax = set_extrema(data, vmin, vmax)
            cbar_ticks = np.linspace(vmin, vmax, nv, endpoint=True)

        if cbar_ax is None:
            cbar_ax = axe

        if cbar_properties is not None:
            cbar = custom_cbar(
                fig, img, cbar_ax, cbar_cax, cbar_ticks, cbar_properties)
        else:
            cbar = fig.colorbar(
                img, ax=cbar_ax, cax=cbar_cax, ticks=cbar_ticks)

        cbar.set_label(data_name)

    # Set frame default limits
    xlim, ylim = set_frame_limits(mesh)
    axe.set_xlim(xlim[0], xlim[1])
    axe.set_ylim(ylim[0], ylim[1])

    # Adding labels
    axe.set_xlabel(x_label)
    axe.set_ylabel(y_label)


def plot2d_scalar_filled_contour(fig, axe, mesh, data,
                                 x_label='', y_label='', data_name='data',
                                 vmin=None, vmax=None, nv=11, levels=None,
                                 colorbar=True, cbar_properties=None,
                                 cbar_ticks=None,
                                 cbar_ax=None, cbar_cax=None,
                                 cmap_name='jet', **kwargs):
    """
    Plot a 2d filled contour of scalar data on triangulation
    (matplotlib tricontourf)

    @param fig (matplotlib.figure) matplotlib figure structure
    @param axe (matplotlib.axes) matplotlib axes on which to draw
    @param mesh (matplotlib.tri.Triangulation or tuple)
    triangular mesh or regular x,y grid
    @param data (numpy.array) scalar data to plot
    @param x_label (string) Name of the x_label (default '')
    @param y_label (string) Name of the y_label (default '')
    @param data_name (str) name of the data to display
    @param vmin (float) Minimal value of data to plot
    @param vmax (float) Maximal value of data to plot
    @param nv (integer) Number of sample for colorbar range (default 11)
    @param levels (np.array) levels used for contours (default None)
    @param colorbar (bool) show colorbar (default: True)
    @param cbar_ticks (list) list of values where to show color bar ticks
    @param cbar_ax (Axes) Parent axes from which space for a new colorbar
    axes will be stolen. If a list of axes is given they will all be resized
    to make room for the colorbar axes.
    @param cbar_cax (Axes) Axes into which the colorbar will be drawn.
    @param cbar_properties (dict) list additional properties of the colorbar
    @param cmap_name (str) colormap name (default: 'jet')
    @param kwargs (dict) rest of optional arguments given to tricontourf
    """
    # Getting cmap object from name
    cmap = cm.get_cmap(name=cmap_name, lut=None)

    # Setting levels with min, max and nv values if levels is not prescribed
    if vmax is not None or vmin is not None:
        assert levels is None
        vmin, vmax = set_extrema(data, vmin, vmax)
        levels = np.linspace(vmin, vmax, nv)

    # Plotting scalar filled contour
    if isinstance(mesh, mtri.triangulation.Triangulation):
        if levels is not None:
            img = axe.tricontourf(mesh, data, levels, cmap=cmap, **kwargs)
        else:
            img = axe.tricontourf(mesh, data, cmap=cmap, **kwargs)
    elif isinstance(mesh, tuple):
        x, y = mesh
        if levels is not None:
            img = axe.contourf(x, y, data, levels, cmap=cmap, **kwargs)
        else:
            img = axe.contourf(x, y, data, cmap=cmap, **kwargs)
    else:
        raise ValueError("Unknown mesh type")

    # Adding color bar if requested
    if colorbar:
        if levels is not None and cbar_ticks is None:
            cbar_ticks = levels

        if cbar_ax is None:
            cbar_ax = axe

        if cbar_properties is not None:
            cbar = custom_cbar(
                fig, img, cbar_ax, cbar_cax, cbar_ticks, cbar_properties)
        else:
            cbar = fig.colorbar(
                img, ax=cbar_ax, cax=cbar_cax, ticks=cbar_ticks)

        cbar.set_label(data_name)

    # Set frame default limits
    xlim, ylim = set_frame_limits(mesh)
    axe.set_xlim(xlim[0], xlim[1])
    axe.set_ylim(ylim[0], ylim[1])

    # Adding labels
    axe.set_xlabel(x_label)
    axe.set_ylabel(y_label)

    return img


def plot2d_scalar_contour(fig, axe, mesh, data,
                          x_label='', y_label='', data_name='data',
                          vmin=None, vmax=None, nv=11, levels=None,
                          colorbar=True, cbar_properties=None, cbar_ticks=None,
                          cbar_ax=None, cbar_cax=None,
                          colors=None, cmap_name='jet', **kwargs):
    """
    Plot a 2d contour of scalar data on triangulation
    (matplotlib tricontour)

    @param fig (matplotlib.figure) matplotlib figure structure
    @param axe (matplotlib.axes) matplotlib axes on which to draw
    @param mesh (matplotlib.tri.Triangulation or tuple)
    triangular mesh or regular x,y grid
    @param data (numpy.array) scalar data to plot
    @param x_label (string) Name of the x_label (default '')
    @param y_label (string) Name of the y_label (default '')
    @param data_name (str) name of the data to display
    @param vmin (float) Minimal value of data to plot
    @param vmax (float) Maximal value of data to plot
    @param nv (integer) Number of sample for colorbar range (default 11)
    @param levels (np.array) levels used for contours (default None)
    @param colorbar (bool) show colorbar (default: True)
    @param cbar_ticks (list) list of values where to show color bar ticks
    @param cbar_ax (Axes) Parent axes from which space for a new colorbar
    axes will be stolen. If a list of axes is given they will all be resized
    to make room for the colorbar axes.
    @param cbar_cax (Axes) Axes into which the colorbar will be drawn.
    @param cbar_properties (dict) list additional properties of the colorbar
    @param colors (str) colors of the contour lines (default: None)
    @param cmap_name (str) colormap of the contour lines (default: 'jet')
    @param kwargs (dict) rest of optional arguments given to tricontour
    """
    # Getting cmap object from name
    cmap = cm.get_cmap(name=cmap_name, lut=None)

    # Setting levels with min, max and nv values if levels is not prescribed
    if vmax is not None or vmin is not None or nv != 11:
        assert levels is None
        vmin, vmax = set_extrema(data, vmin, vmax)
        levels = np.linspace(vmin, vmax, nv)

    # Plotting scalar contour
    if isinstance(mesh, mtri.triangulation.Triangulation):
        if colors is None:
            if levels is not None:
                img = axe.tricontour(mesh, data, levels, cmap=cmap, **kwargs)
            else:
                img = axe.tricontour(mesh, data, cmap=cmap, **kwargs)
        else:
            if levels is not None:
                img = axe.tricontour(mesh, data, levels, colors=colors,
                                     **kwargs)
            else:
                img = axe.tricontour(mesh, data, colors=colors, **kwargs)
    elif isinstance(mesh, tuple):
        x, y = mesh
        if colors is None:
            if levels is not None:
                img = axe.contour(x, y, data, levels, cmap=cmap, **kwargs)
            else:
                img = axe.contour(x, y, data, cmap=cmap, **kwargs)
        else:
            if levels is not None:
                img = axe.contour(x, y, data, levels, colors=colors, **kwargs)
            else:
                img = axe.contour(x, y, data, colors=colors, **kwargs)
    else:
        raise ValueError("Unknown mesh type")

    # Adding color bar if requested
    if colorbar:
        if levels is not None and cbar_ticks is None:
            cbar_ticks = levels

        if cbar_ax is None:
            cbar_ax = axe

        if cbar_properties is not None:
            cbar = custom_cbar(
                fig, img, cbar_ax, cbar_cax, cbar_ticks, cbar_properties)
        else:
            cbar = fig.colorbar(
                img, ax=cbar_ax, cax=cbar_cax, ticks=cbar_ticks)

        cbar.set_label(data_name)

    # Set frame default limits
    xlim, ylim = set_frame_limits(mesh)
    axe.set_xlim(xlim[0], xlim[1])
    axe.set_ylim(ylim[0], ylim[1])

    # Adding labels
    axe.set_xlabel(x_label)
    axe.set_ylabel(y_label)


def plot2d_vectors(fig, axe, mesh, data_x, data_y,
                   x_label='', y_label='', data_name='data',
                   normalize=True, scale=50,
                   vmin=None, vmax=None, nv=None,
                   colorbar=True, cbar_properties=None, cbar_ticks=None,
                   cbar_ax=None, cbar_cax=None,
                   grid_xlim=None, grid_ylim=None, grid_resolution=None,
                   color=None, cmap_name='jet', **kwargs):
    """
    Plot 2d vectors (data_x, data_y) on triangulation
    (matplotlib quiver)

    @param fig (matplotlib.figure) matplotlib figure structure
    @param axe (matplotlib.axes) matplotlib axes on which to draw
    @param mesh (matplotlib.tri.Triangulation or tuple)
    triangular mesh or regular x,y grid
    @param data_x (numpy.array) x component of the vector to plot
    @param data_y (numpy.array) y component of the vector to plot
    @param x_label (string) Name of the x_label (default '')
    @param y_label (string) Name of the y_label (default '')
    @param data_name (str) Name of the data to display
    @param normalize (bool) normalize vectors (default: True)
    @param scale (float) scale of the vectors (default: 50)
    @param vmin (float) Minimal value of data to plot
    @param vmax (float) Maximal value of data to plot
    @param nv (integer) Number of sample for colorbar range (default 10)
    @param colorbar (bool) show colorbar (default: True)
    @param cbar_ticks (list) list of values where to show color bar ticks
    @param cbar_ax (Axes) Parent axes from which space for a new colorbar
    axes will be stolen. If a list of axes is given they will all be resized
    to make room for the colorbar axes.
    @param cbar_cax (Axes) Axes into which the colorbar will be drawn.
    @param cbar_properties (dict) list additional properties of the colorbar
    @param grid_xlim (list(float)) grid limit in x
    @param grid_ylim (list(float)) grid limit in y
    @param grid_resolution (list(float)) number of grid point on x and y
    @param color (str) color of the vectors (default: None)
    @param cmap_name (str) colormap of the vectors (default: 'jet')
    @param kwargs (dict) rest of optional arguments given to quiver
    """
    # Getting cmap object from name
    cmap = cm.get_cmap(name=cmap_name, lut=None)

    # Defining vector locations (on regular grid or triangulation)
    if isinstance(mesh, mtri.triangulation.Triangulation):
        if grid_resolution is None:
            # Triangulation points
            x_i, y_i = mesh.x, mesh.y
        else:
            # Interpolate data_x and data_y on grid
            data_iterpolated, (x_i, y_i) = interpolate_on_grid(
                mesh, [data_x, data_y],
                grid_xlim=grid_xlim, grid_ylim=grid_ylim,
                grid_resolution=grid_resolution)
            data_x = data_iterpolated[0]
            data_y = data_iterpolated[1]

    elif isinstance(mesh, tuple):
        # User provided grid
        x_i, y_i = mesh

    else:
        raise ValueError("Unknown mesh type")

    # Normalize vector
    data_norm = np.sqrt(data_x**2 + data_y**2)
    if normalize:
        data_x = data_x / (data_norm + 1.e-8)
        data_y = data_y / (data_norm + 1.e-8)

    # Plotting colored vectors with cmap
    if color is None:
        img = axe.quiver(x_i, y_i,
                         data_x, data_y,
                         data_norm[:], scale=scale, cmap=cmap,
                         **kwargs)

        # Adding color bar if requested
        if colorbar:
            if vmin is not None or vmax is not None or nv is not None:
                vmin, vmax = set_extrema(data_norm, vmin, vmax)
                img.set_clim(vmin, vmax)

                if nv is not None and cbar_ticks is None:
                    cbar_ticks = np.linspace(vmin, vmax, nv, endpoint=True)

            if cbar_ax is None:
                cbar_ax = axe

            if cbar_properties is not None:
                cbar = custom_cbar(
                    fig, img, cbar_ax, cbar_cax, cbar_ticks, cbar_properties)
            else:
                cbar = fig.colorbar(
                    img, ax=cbar_ax, cax=cbar_cax, ticks=cbar_ticks)

            cbar.set_label(data_name)

    # Plotting unicolored vectors
    else:
        img = axe.quiver(x_i, y_i,
                         data_x, data_y,
                         scale=scale, color=color, **kwargs)

    # Set frame default limits
    xlim, ylim = set_frame_limits(mesh)
    axe.set_xlim(xlim[0], xlim[1])
    axe.set_ylim(ylim[0], ylim[1])

    # Adding labels
    axe.set_xlabel(x_label)
    axe.set_ylabel(y_label)


def plot2d_streamlines(fig, axe, mesh, data_x, data_y,
                       x_label='', y_label='', data_name='data',
                       density=1.5, arrowsize=1.5,
                       vmin=None, vmax=None, nv=None,
                       colorbar=True, cbar_properties=None, cbar_ticks=None,
                       cbar_ax=None, cbar_cax=None,
                       grid_xlim=None, grid_ylim=None,
                       grid_resolution=[500, 500],
                       color=None, cmap_name='jet', linewidth=0.5, **kwargs):
    """
    Plot 2d streamlines based on vectorial data (data_x, data_y)
    (matplotlib streamplot on cartesian subgrid)

    @param fig (matplotlib.figure) matplotlib figure structure
    @param axe (matplotlib.axes) matplotlib axes on which to draw
    @param mesh (matplotlib.tri.Triangulation or tuple)
    triangular mesh or regular x,y grid
    @param data_x (numpy.array) x component of the vector to plot
    @param data_y (numpy.array) y component of the vector to plot
    @param x_label (string) Name of the x_label (default '')
    @param y_label (string) Name of the y_label (default '')
    @param data_name (str) Name of the data to display
    @param density (float) density of the streamlines (default: 1.5)
    @param arrowsize (float) arrow size on strealines (default: 1.5)
    @param vmin (float) Minimal value of data to plot
    @param vmax (float) Maximal value of data to plot
    @param nv (integer) Number of sample for colorbar range (default 10)
    @param colorbar (bool) show colorbar (default: True)
    @param cbar_ticks (list) list of values where to show color bar ticks
    @param cbar_ax (Axes) Parent axes from which space for a new colorbar
    axes will be stolen. If a list of axes is given they will all be resized
    to make room for the colorbar axes.
    @param cbar_cax (Axes) Axes into which the colorbar will be drawn.
    @param cbar_properties (dict) list additional properties of the colorbar
    @param grid_xlim (list(float)) grid limit in x
    @param grid_ylim (list(float)) grid limit in y
    @param grid_resolution (list(float)) number of grid point on x and y
    @param color (str) color of the streamlines (default: None)
    @param cmap_name (str) colormap of the streamlines (default: 'jet')
    @param linewidth (float) linewidth of the streamlines (default: 0.5)
    @param kwargs (dict) rest of optional arguments given to streamplot
    """
    # Getting cmap object from name
    cmap = cm.get_cmap(name=cmap_name, lut=None)

    # Interpolate data_x and data_y on grid
    if isinstance(mesh, mtri.triangulation.Triangulation):
        data, (x_i, y_i) = interpolate_on_grid(
            mesh, [data_x, data_y],
            grid_xlim=grid_xlim, grid_ylim=grid_ylim,
            grid_resolution=grid_resolution)
    elif isinstance(mesh, tuple):
        # User provided grid
        x_i, y_i = mesh
        data = [data_x, data_y]
    else:
        raise ValueError("Unknown mesh type")

    # Compute norm
    data_norm = np.sqrt(data[0]**2 + data[1]**2)

    # Plotting colored streamlines with cmap
    if color is None:
        img = axe.streamplot(x_i, y_i,
                             data[0],
                             data[1],
                             density=(density, density),
                             arrowsize=arrowsize,
                             color=data_norm, cmap=cmap,
                             linewidth=linewidth, **kwargs)

        # Adding color bar if requested
        if colorbar:
            if vmin is not None or vmax is not None or nv is not None:
                vmin, vmax = set_extrema(data_norm, vmin, vmax)
                img.lines.set_clim(vmin, vmax)

                if nv is not None and cbar_ticks is None:
                    cbar_ticks = np.linspace(vmin, vmax, nv, endpoint=True)

            if cbar_ax is None:
                cbar_ax = axe

            if cbar_properties is not None:
                cbar = custom_cbar(
                    fig, img.lines, cbar_ax, cbar_cax, cbar_ticks,
                    cbar_properties)
            else:
                cbar = fig.colorbar(
                    img.lines, ax=cbar_ax, cax=cbar_cax, ticks=cbar_ticks)

            cbar.set_label(data_name)

    # Plotting unicolored streamlines
    else:
        img = axe.streamplot(x_i, y_i,
                             data[0],
                             data[1],
                             density=(density, density),
                             arrowsize=arrowsize,
                             color=color, cmap=cmap,
                             linewidth=linewidth, **kwargs)

    # Set frame default limits
    xlim, ylim = set_frame_limits(mesh)
    axe.set_xlim(xlim[0], xlim[1])
    axe.set_ylim(ylim[0], ylim[1])

    # Adding labels
    axe.set_xlabel(x_label)
    axe.set_ylabel(y_label)


def custom_cbar(fig, img, cbar_ax,
                cbar_cax, cbar_ticks, cbar_properties):
    """
    Adding custom cbar to a figure

    @param fig (matplotlib.figure) matplotlib figure structure
    @param img (matplotlib object) matplotlib mappable object
    @param cbar_ax (Axes) Parent axes from which space for a new colorbar
    axes will be stolen. If a list of axes is given they will all be resized
    to make room for the colorbar axes.
    @param cbar_cax (Axes) Axes into which the colorbar will be drawn.
    @param cbar_properties (dict) list additional properties of the colorbar
    @param cbar_ticks (list) list of values where to show color bar ticks
    (will overwrite ticks contained in cbar_properties)

    @return colorbar (matplotlib.pyplot.colorbar)
    """
    # update default cbar deco with cbar properties
    deco_cbar.update(cbar_properties)

    # if no ticks in cbar_properties overide with levels ticks
    if deco_cbar['ticks'] is None:
        tks = cbar_ticks
    else:
        tks = deco_cbar['ticks']

    # return custop colorbar
    return fig.colorbar(
        img, ax=cbar_ax, cax=cbar_cax,
        orientation=deco_cbar['orientation'],
        fraction=deco_cbar['fraction'],
        pad=deco_cbar['pad'],
        shrink=deco_cbar['shrink'],
        aspect=deco_cbar['aspect'],
        anchor=deco_cbar['anchor'],
        panchor=deco_cbar['panchor'],
        extend=deco_cbar['extend'],
        extendfrac=deco_cbar['extendfrac'],
        extendrect=deco_cbar['extendrect'],
        spacing=deco_cbar['spacing'],
        ticks=tks,
        format=deco_cbar['format'],
        drawedges=deco_cbar['drawedges'])


def plot2d_quadrangle_mesh(fig, axe, x, y, ikle, **kwargs):
    """
    Plotting a 2d spectrum data

    @param fig (matplotlib.figure) matplotlib figure structure
    @param axe (matplotlib.axes) matplotlib axes on which to draw
    @param ikle (numpy.array) Connectivity of quadrangles
    @param x (numpy.array) x coordinates of the points
    @param y (numpy.array) y coordinates of the points
    @param kwargs (dict) rest of optional arguments given to PolyCollection
    """

    # Creating an array of [x[i], y[i]] array
    x_y = np.c_[x, y]
    # Creating list of vertices
    verts = x_y[ikle]
    # Create collection of vertices for matplotlib
    p_c = mpl.collections.PolyCollection(verts, **kwargs)
    # Adding that to the axe
    axe.add_collection(p_c)


def plot2d_spectrum(fig, axe, x, y, ikle, data, **kwargs):
    """
    Plotting a 2d spectrum data

    @param fig (matplotlib.figure) matplotlib figure structure
    @param axe (matplotlib.axes) matplotlib axes on which to draw
    @param ikle (numpy.array) Connectivity of quadrangles
    @param x (numpy.array) x coordinates of the points
    @param y (numpy.array) y coordinates of the points
    @param data (numpy.array) Spectrum data
    @param kwargs (dict) rest of optional arguments given to scalar_map
    """

    # @todo add plot2d_scalar_map arguments
    # Building bogus triangulation by splitting quadrangles in two
    nelem, _ = ikle.shape
    triangles = np.zeros((nelem*2, 3))
    # Splitting quadrangles
    for i, elem in enumerate(ikle):
        i_1, i_2, i_3, i_4 = elem
        triangles[i] = [i_1, i_2, i_3]
        triangles[i+nelem] = [i_1, i_3, i_4]

    tri = mtri.Triangulation(x, y, triangles)

    # Potting our result
    plot2d_scalar_map(fig, axe, tri, data, **kwargs)
