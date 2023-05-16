"""
Contains the class Selafins
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import numpy as np

# ~~> dependencies towards other pytel/modules
from data_manip.formats.selafin import Selafin
from utils.progressbar import ProgressBar
from utils.exceptions import TelemacException


class Selafins(object):
    """
    Class to handle multiple Serafin files
    It is to handle merge of files
    """

    def __init__(self):
        """
        Initialisation function
        """
        self.slfs = []
        self.slf = None
        self.suite = True
        self.merge = True

    def add(self, file_name):
        """
        Adding a Selafin file to the class

        @param file_name (string) Name of the Serafin file
        """
        slf = Selafin(file_name)
        if self.slf is None:
            self.slf = slf
        self.slfs.append(slf)
        self.suite = self.is_suite()
        self.merge = self.is_merge()

    def is_suite(self):
        """
        Returns true if all the files have the same variables as the first one

        @return (boolean) value of the test
        """
        same = True
        for slf in self.slfs[1:]:
            for var in slf.varnames:
                same = same and (var in self.slfs[0].varnames)
            for var in slf.cldnames:
                same = same and (var in self.slfs[0].cldnames)
        return same

    def is_merge(self):
        """
        Returns true if all the files have the same time steps

        @return (boolean) value of the test
        """
        same = True
        for slf in self.slfs[1:]:
            mmax = 1.e-5 + np.max(slf.tags['times']) +\
                     np.max(self.slfs[0].tags['times'])
            # /!\ max always positive
            accuracy = np.power(10.0, -5+np.floor(np.log10(mmax)))
            if len(slf.tags['times']) != len(self.slfs[0].tags['times']):
                same = False
            else:
                for time1, time2 in zip(slf.tags['times'],
                                        self.slfs[0].tags['times']):
                    same = same and (accuracy > (time1-time2))
        return same

    def pop(self, index=0):
        """
        Returns and remove a file from the list of files.
        If the index is given the file at that index is taken.

        @param index (int) index of the file to take (default=0)

        @return (Selafin) A selafin object
        """
        index = max(0, min(index, len(self.slfs)-1))
        return self.slfs.pop(self.slfs[index])

    # TODO: files also have to have the same header
    def put_content(self, file_name, showbar=True):
        """
        Write the merge/suite of all the files contained in the class into a
        single file

        @param file_name (string) Name of the output file
        @param showbar (boolean) Display progression bar
        """
        if self.suite and self.merge:
            if len(self.slfs) == 2:  # /!\ difference only between two files
                self.slf.fole.update({'hook': open(file_name, 'wb')})
                self.slf.fole['name'] = file_name
                ibar = 0
                if showbar:
                    pbar = \
                       ProgressBar(maxval=len(self.slf.tags['times'])).start()
                self.slf.append_header_slf()
                for itime in range(len(self.slf.tags['times'])):
                    ibar += 1
                    self.slf.append_core_time_slf(itime)
                    self.slf.append_core_vars_slf(self.slf.get_values(itime)
                                                  - self.slfs[1]
                                                  .get_values(itime))
                    if showbar:
                        pbar.update(ibar)
                if showbar:
                    pbar.finish()
            else:
                self.slf.put_content(file_name)  # just a copy
        elif self.suite:
            self.slf.fole.update({'hook': open(file_name, 'wb')})
            self.slf.fole['name'] = file_name
            ibar = 0
            if showbar:
                pbar = ProgressBar(maxval=len(self.slf.tags['times'])).start()
            self.slf.append_header_slf()
            for itime in range(len(self.slf.tags['times'])):
                ibar += 1
                time = self.slf.tags['times'][itime]
                self.slf.append_core_time_slf(itime)
                self.slf.append_core_vars_slf(self.slf.get_values(itime))
                if showbar:
                    pbar.update(ibar)
            if showbar:
                pbar.finish()
            for slf in self.slfs:
                slf.fole = self.slf.fole
                ibar = 0
                if showbar:
                    pbar = ProgressBar(maxval=len(slf.tags['times'])).start()
                for itime in range(len(slf.tags['times'])):
                    if slf.tags['times'][itime] > time:
                        ibar += 1
                        time = slf.tags['times'][itime]
                        slf.append_core_time_slf(itime)
                        slf.append_core_vars_slf(slf.get_values(itime))
                        if showbar:
                            pbar.update(ibar)
                if showbar:
                    pbar.finish()
            self.slf.fole['hook'].close()
        elif self.merge:
            self.slf.fole.update({'hook': open(file_name, 'wb')})
            self.slf.fole['name'] = file_name
            for slf in self.slfs[1:]:
                slf.fole = self.slf.fole
                idvars = []
                for var in range(len(slf.varnames)):
                    if var not in self.slf.varnames:
                        idvars.append(var)
                        self.slf.varnames.append(slf.varnames[var])
                        self.slf.varunits.append(slf.varunits[var])
                for var in range(len(slf.cldnames)):
                    if var not in self.slf.cldnames:
                        idvars.append(var+slf.nbv1)
                        self.slf.cldnames.append(slf.cldnames[var])
                        self.slf.cldunits.append(slf.cldunits[var])
                slf.varindex = idvars
                self.slf.nbv1 = len(self.slf.varnames)
                self.slf.nbv2 = len(self.slf.cldnames)
            ibar = 0
            if showbar:
                pbar = ProgressBar(maxval=len(self.slf.tags['times'])).start()
            self.slf.append_header_slf()
            for itime in range(len(self.slf.tags['times'])):
                ibar += 1
                self.slf.append_core_time_slf(itime)
                self.slf.append_core_vars_slf(self.slf.get_values(itime))
                for slf in self.slfs[1:]:
                    self.slf.append_core_vars_slf(slf.get_values(itime))
                if showbar:
                    pbar.update(ibar)
            if showbar:
                pbar.finish()
            self.slf.fole['hook'].close()
        elif len(self.slf.tags['times']) == 1 and len(self.slfs) == 2:
            # self.slf will be distributed over the time frames
            # of the scond other
            self.slf.fole.update({'hook': open(file_name, 'wb')})
            self.slf.fole['name'] = file_name
            slf = self.slfs[1]
            slf.fole = self.slf.fole
            idvars = []
            for var in range(len(slf.varnames)):
                if var not in self.slf.varnames:
                    idvars.append(var)
                    self.slf.varnames.append(slf.varnames[var])
                    self.slf.varunits.append(slf.varunits[var])
            for var in range(len(slf.cldnames)):
                if var not in self.slf.cldnames:
                    idvars.append(var+slf.nbv1)
                    self.slf.cldnames.append(slf.cldnames[var])
                    self.slf.cldunits.append(slf.cldunits[var])
            slf.varindex = idvars
            self.slf.nbv1 = len(self.slf.varnames)
            self.slf.nbv2 = len(self.slf.cldnames)
            ibar = 0
            if showbar:
                pbar = ProgressBar(maxval=len(slf.tags['times'])).start()
            self.slf.append_header_slf()
            for itime in range(len(slf.tags['times'])):
                ibar += 1
                slf.append_core_time_slf(slf, itime)
                self.slf.append_core_vars_slf(self.slf.get_values(0))
                slf.append_core_vars_slf(slf.get_values(itime))
                if showbar:
                    pbar.update(ibar)
            if showbar:
                pbar.finish()
            self.slf.fole['hook'].close()
        else:
            raise TelemacException(
                 "Does not know how to merge your files. Try either:\n"
                 "    + to make sure your files have the same time support"
                 "    + to make sure your files have the same variables")

    def __del__(self):
        """
        Deleting the object
        """
        if self.slf is not None:
            del self.slf
        if self.slfs != []:
            for slf in self.slfs:
                del slf
