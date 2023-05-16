"""
Contains the class PARAFINS
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
from os import path, getcwd
from copy import deepcopy
import glob
# ~~> dependencies towards other pytel/modules
from data_manip.formats.selafins import Selafins
from data_manip.formats.selafin import Selafin
from utils.progressbar import ProgressBar
from utils.exceptions import TelemacException

import numpy as np


class Parafins(Selafins):
    """
    Class to handle partitionned serafin files
    """

    def __init__(self, file_name, root=None):
        """
        Initialise the class

        @param file_name (string) Name of the global selafin file
        @param root (string) Name of the non partitionned file

        """
        Selafins.__init__(self)
        # ~~> The main slf is based on the header of the GEO file
        self.slf = Selafin(file_name)
        self.file = self.slf.file
        self.fole = self.slf.fole
        if root is not None:
            # ~~> Loading the individual headers
            self.add_root(root)
            # ~~> Making sure there are all inter-compatible
            if self.suite and self.merge:
                self.slf.tags = self.slfs[0].tags
                self.slf.nbv1 = self.slfs[0].nbv1
                self.slf.varnames = self.slfs[0].varnames
                self.slf.varunits = self.slfs[0].varunits
                self.slf.nbv2 = self.slfs[0].nbv2
                self.slf.cldnames = self.slfs[0].cldnames
                self.slf.cldunits = self.slfs[0].cldunits
                self.slf.nvar = self.slf.nbv1 + self.slf.nbv2
                self.slf.varindex = range(self.slf.nvar)
            else:
                raise TelemacException(
                        "... Incompatibilities between files "
                        "for {}".format(path.basename(root)))

            # ~~> Create the corresponding map
            self.map_poin = np.zeros(self.slf.npoin3, dtype=np.int)
            for i, slf in zip(range(len(self.slfs)), self.slfs):
                self.map_poin[slf.ipob3-1] = i

    def alter_endian(self):
        """
        Alter Endianess
        """
        self.slf.alter_endian()

    def alter_float(self):
        """
        Alter float precision
        """
        self.slf.alter_float()

    def add_root(self, root):
        """
        Add a list of partitioned files

        @param root (string) Name of the non partitionned file
        """
        # ~~> list all entries
        diroot = path.dirname(root)
        if path.dirname(root).strip() == '':
            diroot = getcwd()
        root = path.join(diroot, path.basename(root))
        slfnames = glob.glob(root+'?????-?????')
        # ~~> match expression
        if slfnames == []:
            print("... Could not find any sub-files to the root: "+root)
            return []
        npsize = len(slfnames)
        for nptime in range(npsize):
            fle = root+'{0:05d}-{1:05d}'.format(npsize-1, nptime)
            if fle not in slfnames:
                print("... Could not find the following sub-file in "
                      "the list: "+fle)
                return []
        print('      +> Reading the header from the following partitions:')
        for fle in sorted(slfnames):
            print('         ~> '+path.basename(fle))
            slf = Selafin(fle)
            self.slfs.append(slf)
        return

    def get_palues(self, time):
        """
        get values for a given time step on the global mesh

        @param time (int) Time step

        @return (np.array) Containing value for each variable
        """
        _, fsize = self.file['float']
        if len(self.slfs) > 0:
            if fsize == 4:
                varsor = np.zeros((self.slf.nbv1+self.slf.nbv2,
                                   self.slf.npoin3),
                                  dtype=np.float32)
            else:
                varsor = np.zeros((self.slf.nbv1+self.slf.nbv2,
                                   self.slf.npoin3),
                                  dtype=np.float64)
            for slf in self.slfs:
                varsub = slf.get_values(time)
                for ivar in range(self.slf.nvar):
                    varsor[ivar][slf.ipob3-1] = varsub[ivar]
        else:
            varsor = self.slf.get_values(time)
        return varsor

    def get_series(self, nodes, vars_indexes=None):
        """
        get series (value on a point for all time step) for the global mesh

        @param nodes (list) List of nodes to extract series for
        @param vars_indexes (list) List of variables for which to extract data

        @return (np.array) containing the series
        """
        _, fsize = self.file['float']
        if vars_indexes is None:
            vars_indexes = self.slf.varindex
        if len(self.slfs) > 0:
            if fsize == 4:
                z = np.zeros((len(vars_indexes), len(nodes),
                              len(self.slf.tags['cores'])),
                             dtype=np.float32)
            else:
                z = np.zeros((len(vars_indexes), len(nodes),
                              len(self.slf.tags['cores'])),
                             dtype=np.float64)
            mproc = self.map_poin[nodes]
            print('      +> Extracting time series from the following'
                  ' partitions:')
            for islf, slf in zip(range(len(self.slfs)), self.slfs):
                if islf not in mproc:
                    continue
                # ~~> Filter the list of nodes according to sub ipobo
                sub_g_nodes = np.compress(mproc == islf,
                                          np.array(list(zip(range(len(nodes)),
                                                            nodes)),
                                                   dtype=[('r', int),
                                                          ('n', int)]))
                # /!\ why does this work in a sorted way ?
                sub_l_nodes = np.searchsorted(np.sort(slf.ipob2),
                                              sub_g_nodes['n']) + 1
                print('         ~> '+str(len(sub_l_nodes)) +
                      ' nodes from partition '+str(islf))
                # ~~> Get the series from individual sub-files
                subz = slf.get_series(sub_l_nodes, vars_indexes)
                # ~~> Reorder according to original list of nodes
                for ivar in range(len(vars_indexes)):
                    z[ivar][sub_g_nodes['r']] = subz[ivar]
            return z
        else:
            return self.slf.get_series(nodes)  # ,vars_indexes not necessary

    # TODO: files also have to have the same header
    def put_content(self, file_name, showbar=True):
        """
        Write global content into file_name

        @param file_name (string) Name of the output file
        @param showbar (boolean) display a showbar (default=True)
        """
        self.slf.fole.update({'hook': open(file_name, 'wb')})
        ibar = 0
        if showbar:
            pbar = ProgressBar(maxval=len(self.slf.tags['times'])).start()
        self.slf.append_header_slf()
        # ~~> Time stepping
        for itime in range(len(self.slf.tags['times'])):
            ibar += 1
            self.slf.append_core_time_slf(itime)  # Time stamps
            self.slf.append_core_vars_slf(self.get_palues(itime))
            if showbar:
                pbar.update(ibar)
        self.slf.fole['hook'].close()
        if showbar:
            pbar.finish()

    def cut_content(self, root):
        """
        Apply existing partition to other files

        @param root (string) Name of the global geometry file
        """
        islf = Selafin(root)
        print('      +> Writing the core of the following partitions:')
        # TODO: do this loop in python parallel with a
        # bottle neck at islf.get_values(t)
        for slf in self.slfs:
            sub = deepcopy(slf)
            sub.fole.update({'endian': islf.fole['endian']})
            sub.fole.update({'float': islf.fole['float']})
            _, fsize = islf.fole['float']
            # ~~ Conversion to local islf ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # you know the geom is 2D
            # keep the ipobo2, X2D, Y2D and ikle2 and replace the rest
            sub.nplan = islf.nplan
            sub.iparam = islf.iparam

            # ~~> matching partitions
            if self.slf.npoin2 == islf.npoin2:
                sub.ndp2 = islf.ndp2
                sub.ndp3 = islf.ndp3
                sub.npoin3 = sub.npoin2*islf.nplan
                if islf.nplan > 1:
                    sub.nelem3 = sub.nelem2*(islf.nplan-1)
                    tmp = np.repeat(sub.ipob2, islf.nplan)\
                            .reshape((sub.npoin2, islf.nplan))
                    sub.ipob3 = np.ravel(np.add(tmp,
                                         sub.npoin2*np.arange(islf.nplan)).t)
                    sub.ikle3 = \
                        np.repeat(sub.npoin2*np.arange(islf.nplan-1),
                                  sub.nelem2*islf.ndp3)\
                          .reshape((sub.nelem2*(islf.nplan-1), islf.ndp3)) + \
                        np.tile(np.add(np.tile(sub.ikle2, 2),
                                       np.repeat(sub.npoin2*np.arange(2),
                                                 sub.ndp2)),
                                (islf.nplan-1, 1))
                else:
                    sub.nelem3 = sub.nelem2
                    sub.ipob3 = sub.ipob2
                    sub.ikle3 = sub.ikle2
                indices = sub.ipob2-1

            # ~~> filtered partitions / non reversable !
            else:
                # Intersection with the local domain (sub of slfs):
                # the intersection could be empty
                intsect = np.in1d(islf.ipob2, np.sort(sub.ipob2))
                # The new compressed ipobo
                sub.ipob2 = np.compress(intsect, islf.ipob2)
                # Those node numbers of the islf you keep ~> converting
                # True/False into indices
                indices = np.arange(len(islf.ipob2), dtype=np.int)[intsect]
                # Set the array that only includes elements of islf.ikle2
                # with at least two nodes in the subdomain
                ikles2_1d = np.in1d(islf.ikle2, np.sort(indices))\
                              .reshape(islf.nelem2, islf.ndp2)
                gkle2 = islf.ikle2[np.where(np.sum(ikles2_1d, axis=1) ==
                                            islf.ndp2)]
                # re-numbering ikle2 as a local connectivity matrix
                knolg = np.unique(np.ravel(gkle2))
                knogl = dict(zip(knolg, range(len(knolg))))
                sub.ikle2 = - np.ones_like(gkle2, dtype=np.int)
                for k in range(len(gkle2)):
                    for k_i in range(islf.ndp2):
                        # /!\ sub.ikle2 has a local numbering,
                        # fit to the boundary elements
                        sub.ikle2[k][k_i] = knogl[gkle2[k][k_i]]
                # Set the remaining integers
                sub.npoin2 = len(sub.ipob2)
                sub.nelem2 = len(sub.ikle2)
                sub.ndp2 = islf.ndp2
                sub.ndp3 = islf.ndp3
                sub.npoin3 = sub.npoin2*islf.nplan
                if islf.nplan > 1:
                    sub.nelem3 = sub.nelem2*(islf.nplan-1)
                    tmp = np.repeat(sub.ipob2, islf.nplan)\
                            .reshape((sub.npoin2, islf.nplan))
                    sub.ipob3 = np.ravel(np.add(tmp,
                                         sub.npoin2*np.arange(islf.nplan)).T)
                    sub.ikle3 = \
                        np.repeat(sub.npoin2*np.arange(islf.nplan-1),
                                  sub.nelem2*islf.ndp3)\
                          .reshape((sub.nelem2*(islf.nplan-1), islf.ndp3)) + \
                        np.tile(np.add(np.tile(sub.ikle2, 2),
                                       np.repeat(sub.npoin2*np.arange(2),
                                                 sub.ndp2)),
                                (islf.nplan-1, 1))
                else:
                    sub.nelem3 = sub.nelem2
                    sub.ipob3 = sub.ipob2
                    sub.ikle3 = sub.ikle2
                # Set the remaining floats
                sub.meshx = islf.meshx[indices]
                sub.meshy = islf.meshy[indices]

            # ~~ Addition of variables from islf ~~~~~~~~~~~~~~~~~~~~~~
            sub.title = islf.title
            sub.datetime = islf.datetime
            sub.varnames = islf.varnames
            sub.varunits = islf.varunits
            sub.cldnames = islf.cldnames
            sub.cldunits = islf.cldunits
            sub.nbv1 = islf.nbv1
            sub.nbv2 = islf.nbv2
            sub.nvar = sub.nbv1+sub.nbv2
            # ~~ put_content ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            file_name = path.basename(sub.file['name'])\
                .replace(self.slf.file['name'], root)
            sub.fole.update({'hook': open(file_name, 'wb')})
            sub.append_header_slf()
            print('         ~> '+path.basename(file_name))
            # ~~> Time stepping
            sub.tags['times'] = islf.tags['times']
            if fsize == 4:
                varsors = np.zeros((sub.nvar, sub.npoin3), dtype=np.float32)
            else:
                varsors = np.zeros((sub.nvar, sub.npoin3), dtype=np.float64)
            for itime in range(len(islf.tags['times'])):
                sub.append_core_time_slf(itime)  # Time stamps
                for ivar, var in zip(range(sub.nvar), islf.get_values(itime)):
                    for n in range(sub.nplan):
                        varsors[ivar][n*sub.npoin2:(n+1)*sub.npoin2] = \
                              var[n*islf.npoin2:(n+1)*islf.npoin2][indices]
                sub.append_core_vars_slf(varsors)
            sub.fole['hook'].close()
