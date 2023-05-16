r"""@author TELEMAC-MASCARET Consortium

"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
# ~~> dependencies towards other modules
# ~~> dependencies towards other modules
from pretel.chop_selafin import ChopSelafin

class AlterSelafin(ChopSelafin):

    alter_zp = 0
    alter_zm = 1
    alter_z_names = []

    def __init__(self, f, times=None, vrs=None, root=None):
        ChopSelafin.__init__(self, f, root)
        if vrs != None:
            ChopSelafin.update_vars(self, vrs)
        if times != None:
            tfrom, tstep, tstop = times
            ChopSelafin.update_times(self, tfrom, tstep, tstop)

    def alter_title(self, title=None):
        if title != None:
            self.slf.title = (title + 80*' ')[0:80]

    def alter_datetime(self, date=None, time=None):
        if date != None:
            for i in range(3):
                self.slf.datetime[i] = int(date[2-i])
            self.slf.iparam[9] = 1
        if time != None:
            self.slf.datetime[3] = int(time[0])
            if len(time) > 1:
                self.slf.datetime[4] = int(time[1])
            if len(time) > 2:
                self.slf.datetime[5] = int(time[2])
            self.slf.iparam[9] = 1

    def alter_vars(self, vrs=None):
        if vrs != None:
            for var_n in vrs.split(':'):
                v, n = var_n.split('=')
                for ivar in range(len(self.slf.varnames)):
                    if v.lower() in self.slf.varnames[ivar].lower():
                        self.slf.varnames[ivar] = n.upper()
                for ivar in range(len(self.slf.cldnames)):
                    if v.lower() in self.slf.cldnames[ivar].lower():
                        self.slf.cldnames[ivar] = n.upper()

    def switch_vars(self):
        x = self.slf.varindex[0:self.slf.nbv1]
        self.slf.varindex = self.slf.varindex[self.slf.nbv1:]
        self.slf.varindex.extend(x)
        x = self.slf.nbv1
        self.slf.nbv1 = self.slf.nbv2
        self.slf.nbv2 = x
        x = self.slf.varnames
        self.slf.varnames = self.slf.cldnames
        self.slf.cldnames = x
        x = self.slf.varunits
        self.slf.varunits = self.slf.cldunits
        self.slf.cldunits = x

    def alter_mesh(self, m_x=1, p_x=0, m_y=1, p_y=0):
        self.slf.meshx = m_x * self.slf.meshx + p_x
        self.slf.meshy = m_y * self.slf.meshy + p_y

    def alter_times(self, m_t=1, p_t=0):
        self.slf.tags['times'] = m_t * self.slf.tags['times'] + p_t

    def alter_values(self, vrs=None, m_z=1, p_z=0):
        self.slf.alter_values(vrs, m_z, p_z)
