#!/usr/bin/env python3
"""
Functions to interpolation bathymetry into a file
"""
from functools import wraps
import numpy as np
from scipy import interpolate

def lazy_property(fn):
    """
    Version of lazy_property by John Huang.

    Decorator used to cache property results into dictionary.
    The cache can be cleared using clean_lazy_properties.
    """
    cache_name = '_lazy_properties'
    attr_name = fn.__name__

    def get_cache(instance):
        if not hasattr(instance, cache_name):
            setattr(instance, cache_name, {})
        return getattr(instance, cache_name)

    @property
    @wraps(fn)
    def get_attr(self):
        cache = get_cache(self)
        if attr_name not in cache:
            cache[attr_name] = fn(self)
        return cache[attr_name]

    return get_attr


class Interp2D(object):
    """
    2D interpolator

    @param xyz (numpy.array) x, y coordonates of points and
    corresponding values z = f(x, y) to interpolate
    @param kind (str) Interpolation method. Can be:
        - 'nearest'
        - 'linear'
        - 'cubic'

    """
    def __init__(self, xyz, kind='nearest', **kwargs):

        self.x = np.array(xyz[:, 0])
        self.y = np.array(xyz[:, 1])
        self.z = np.array(xyz[:, 2])
        self.kind = kind
        self.kwargs = kwargs

    @lazy_property
    def _template_interp(self):
        """
        Construct template interpolator function based on kind
        """
        if self.kind == 'nearest':
            template = interpolate.NearestNDInterpolator

        elif self.kind == 'linear':
            template = interpolate.LinearNDInterpolator

        elif self.kind == 'cubic':
            template = interpolate.CloughTocher2DInterpolator

        else:
            raise ValueError("Unknown interpolation method")

        return template

    def interpolate(self, x, y):
        """
        Interpolate points.

        @param x (array of shape (m,)) x-coordinates of points to interpolate
        @param y (array of shape (n,)) y-coordinates of points to interpolate
        @return values (array of shape (m, n)) interpolated values of points.
        """
        points = np.column_stack((self.x, self.y))
        interp = self._template_interp(points, self.z, **self.kwargs)

        return interp(np.column_stack((x, y)))

    def __call__(self, x, y):
        """
        Interpolate in the style of LinearNDInterpolator.

        @param x (array of shape (m,)) x-coordinates of points to interpolate
        @param y (array of shape (n,)) y-coordinates of points to interpolate
        @return values (array of shape (m, n)) interpolated values of points.
        """
        return self.interpolate(x, y)


class Interp2DChunk(Interp2D):
    """
    Reconstruction of interpolation for 2d applications.
    This class is used to avoid any memory errors due to interpolation
    of large numbers of points.

    Built for use for extremely large point clouds. Interpolation
    is partitioned into automatic control parameters px, py, pe, blockpts.
    The scipy implementation of interpolation functions has memory problems
    for large point clouds. This class divides the problem into several
    smaller partitions.

    @param px (int)
        Number of partitions in x-direction. If None, a default is calculated
        according to the number of blockpts
    @param py (int)
        Number of partitions in y-direction. If None, a default is calculated
        according to the number of blockpts.
    @param pe (float)
        Proportion of block length to overlap on other blocks.
        For example, if pe=0.25, the block will be extended 25% on both the
        left and right sides of px to overlap on successive blocks.
    @param blockpts (int)
        Approximate number of interpolation points within each partition block.
        Defaults to 300*300. blockpts is used to automatically size either
        px or py if these are set to None.

    """
    def __init__(self, xyz, kind='nearest',
                 px=None, py=None, pe=0.5, blockpts=300*300,
                 **kwargs):

        Interp2D.__init__(self, xyz, kind=kind, **kwargs)

        self.values = self.z
        self.px = px
        self.py = py
        self.pe = pe
        self.blockpts = blockpts
        self._set_partitions()


    def _set_partitions(self):
        """ Calculate the number of partitions to use in data set"""
        ptnum = len(self.x)
        blockpts = self.blockpts

        blocknum = ptnum / blockpts + 1
        if self.px is None:
            if self.py is None:
                self.px = int(np.sqrt(blocknum))
                self.py = int(blocknum / self.px)
            else:
                self.px = int(blocknum / self.py)

        if self.py is None:
            self.py = int(blocknum / self.px)


        self.px = max(self.px, 1)
        self.py = max(self.py, 1)

        self.xmax = np.max(self.x)
        self.xmin = np.min(self.x)
        self.xlen = self.xmax - self.xmin
        self.xp = self.xlen / self.px       # block x length
        self.xe = self.xp * self.pe         # block x overlap length

        self.ymax = np.max(self.y)
        self.ymin = np.min(self.y)
        self.ylen = self.ymax - self.ymin
        self.yp = self.ylen / self.py       # block y length
        self.ye = self.yp * self.pe         # block y overlap length


        xfudge = (self.xmax - self.xmin) / 1000.
        yfudge = (self.ymax - self.ymin) / 1000.

        # Construct block upper/lower limits
        xl = self.xmin - xfudge
        xu = self.xmax + xfudge
        yl = self.ymin - yfudge
        yu = self.ymax + yfudge

        # Construct blocks
        self.xblocks = np.linspace(xl, xu, self.px + 1)
        self.yblocks = np.linspace(yl, yu, self.py + 1)


    def _choose_block(self, x, y):
        """
        Calculate which interpolation block to use for the given
        coordinates (x, y)

        @returns xindex (int array of shape (N,)) index locations
        for x-dimension of blocks
        @returns yindex (int array of shape (N,)) index locations
        for y-dimension of blocks
        """
        xindex = np.searchsorted(self.xblocks, x) - 1
        yindex = np.searchsorted(self.yblocks, y) - 1
        return xindex, yindex

    @lazy_property
    def _interpolators(self):
        """
        Construct interpolators for every block.

        - 0 dimension corresponds to x data.
        - 1 dimension corresponds to y data.
        """

        # Bounds of block interpolation points
        xl_arr = self.xblocks[0:-1] - self.xe
        xu_arr = self.xblocks[1:]  + self.xe

        yl_arr = self.yblocks[0:-1] - self.ye
        yu_arr = self.yblocks[1:] + self.ye

        compteur = 0

        # Loop through all block boundaries and construct interpolators.
        interpolators = []
        for (xl, xu) in zip(xl_arr, xu_arr):
            interpx = []
            for (yl, yu) in zip(yl_arr, yu_arr):

                #Set original data partition
                ix0 = np.logical_and(xl <= self.x, self.x <= xu)
                iy0 = np.logical_and(yl <= self.y, self.y <= yu)
                index1 = np.logical_and(ix0, iy0)
                x0 = self.x[index1]
                y0 = self.y[index1]
                z0 = self.z[index1]
                points = np.column_stack((x0, y0))
                try:
                    interp1 = self._template_interp(points, z0, **self.kwargs)
                    compteur += 1
                    print("Interpolateur ", compteur)
                    interpx.append(interp1)
                except ValueError:
                    interpx.append(None)
            interpolators.append(interpx)
        return interpolators

    def interpolate(self, x, y):
        """Interpolate points.

        @param x (array of shape (m,)) x-coordinates of points to interpolate
        @param y (array of shape (n,)) y-coordinates of points to interpolate
        @return values (array of shape (m, n)) interpolated values of points.
        """
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        xlen = len(x)

        # Property shape the result
        shape = list(self.z.shape)
        shape[0] = xlen
        result = np.empty(shape)
        result[:] = np.nan

        # Loop through all block boundaries and send points to the block's
        # corresponding interpolator.

        xindex, yindex = self._choose_block(x, y)
        for ix in range(self.px):
            for iy in range(self.py):
                index1 = xindex == ix
                index2 = yindex == iy
                index = np.logical_and(index1, index2)
                interp = self._interpolators[ix][iy]

                points = np.column_stack((x[index], y[index]))
                if len(points) > 0:
                    result[index] = interp(points)
        return result
