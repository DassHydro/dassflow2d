"""
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
import numpy as np
from datetime import datetime
import math
import time

from urllib.response import addinfourl
from urllib.request import HTTPRedirectHandler, build_opener
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
from http.client import BadStatusLine

import traceback
# ~~> dependencies towards other modules
from data_manip.formats.selafin import Selafin
from utils.progressbar import ProgressBar
from utils.exceptions import TelemacException
from scipy.io import netcdf

try:
    import json
except ImportError:
    import simplejson as json

# When using the following user ID to access the ECMWF computing
#    facilities it is with the understanding that you acknowledge
#    the terms and conditions detailed at:
#    http://www.ecmwf.int/en/forecasts/software-and-tools
#
# username: s.bourban@hrwallingford.com
# password: Lxtc14
#
CONFIG = {
     "url": "https://api.ecmwf.int/v1",
     "key": "70f6a4499dddb7d17545f9bd3cf5ef3f",
     "email": "s.bourban@hrwallingford.com",
     "password": "Lxtc14",
         }
# ECMWF Re-Analysis, Keys to retrieve requests
#
#  'dataset' (those if general license), "interim" being the default:
#     era15  - ECMWF Global Reanalysis Data - ERA-15 (Jan 1979 - Dec 1993)
#     era20c  - Reanalysis of the 20th-century using surface observations only
#               (Jan 1900 - Dec 2010)
#     era20cmv0  - ERA-20CM: Ensemble of climate model integrations
#                   (Experimental version)
#     era40  - ECMWF Global Reanalysis Data - ERA-40 (Sep 1957 - Aug 2002)
#     eraclim  - ERA-20CM: Ensemble of climate model integrations
#     icoads  - ICOADS v2.5.1 with interpolated 20CR feedback
#     interim  - ECMWF Global Reanalysis Data -
#           ERA Interim (Jan 1979 - present)
#     ispd  - ISPD v2.2
#     yotc  - YOTC (Year of Tropical Convection)
#  'step', "6" being the default:
#     "24/to/120/by/24", ...
#  'number'  : "all",
#  'levtype' : "sl",
#  'date'    : "20071001",
#  'time'    : "00",
#  'origin'  : "all",
#  'type'    : "pf",
#  'param'   : "tp",
#  'area'    : "70/-130/30/-60",
#  'grid'    : "2/2",
#  'target'  : "data.grib" # wil not be used anymore ...
#              since directly into Selafin
#
# _____                  ___________________________________________
# ____/ API Classes /__________________________________________/
#
# (C) Copyright 2012-2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental
# organisation nor
# does it submit to any jurisdiction.
#


class RetryError(Exception):
    def __init__(self, code, text):
        self.code = code
        self.text = text

    def __str__(self):
        return "%d %s" % (self.code, self.text)


class APIException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def robust(func):

    def wrapped(*args, **kwargs):
        tries = 0
        while True:
            try:
                return func(*args, **kwargs)
            except HTTPError as excpt:
                print("WARNING: httplib2.HTTPError received %s" % (str(excpt)))
                if excpt.code < 500:
                    raise
                tries += 1
                if tries > 10:
                    raise
                time.sleep(60)
            except BadStatusLine as excpt:
                print("WARNING: httplib.BadStatusLine received %s" %
                      (str(excpt)))
                tries += 1
                if tries > 10:
                    raise
                time.sleep(60)
            except URLError as excpt:
                print("WARNING: httplib2.URLError received %s %s" %
                      (str(excpt.errno), str(excpt)))
                tries += 1
                if tries > 10:
                    raise
                time.sleep(60)
            except APIException:
                raise
            except RetryError as excpt:
                print("WARNING: HTTP received %s" % (str(excpt.code)))
                print(excpt.text)
                tries += 1
                if tries > 10:
                    raise
                time.sleep(60)
            except Exception as e:
                print("Unexpected error: %s" % (str(sys.exc_info()[0])))
                print(traceback.format_exc())
                raise

    return wrapped


SAY = True
URL = ""


class Ignore303(HTTPRedirectHandler):

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if code in [301, 302]:
            # We want the posts to work even if we are redirected
            if code == 301:
                # TODO: URL Seems undefined
                global SAY, URL
                if SAY:
                    o = req.get_full_url()
                    n = newurl
                    while o != URL and len(o) and len(n) and o[-1] == n[-1]:
                        o = o[0:-1]
                        n = n[0:-1]
                    print()
                    print("*** ECMWF API has moved")
                    print("***   OLD: %s" % (o))
                    print("***   NEW: %s" % (n))
                    print("*** Please update your ~/.ecmwfapirc file")
                    print()
                    SAY = False
            data = None
            if req.has_data():
                data = req.get_data()
            return Request(newurl, data=data, headers=req.headers,
                           origin_req_host=req.get_origin_req_host(),
                           unverifiable=True)
        return None

    def http_error_303(self, req, fp, code, msg, headers):
        infourl = addinfourl(fp, headers, req.get_full_url())
        infourl.status = code
        infourl.code = code
        return infourl


class Connection(object):

    def __init__(self, email=None, key=None, verbose=False, quiet=False):
        self.email = email
        self.key = key
        self.retry = 5
        self.location = None
        self.done = False
        self.value = True
        self.offset = 0
        self.verbose = verbose
        self.quiet = quiet
        self.status = None
        self.last = None

    @robust
    def call(self, url, payload=None, method="GET"):

        if self.verbose:
            print(method + ' ' + url)

        headers = {"Accept": "application/json",
                   "From": self.email,
                   "X-ECMWF-KEY": self.key}

        opener = build_opener(Ignore303)

        data = None
        if payload is not None:
            data = json.dumps(payload)
            data = data.encode('utf-8')
            headers["Content-Type"] = "application/json"

        url = "%s?offset=%d&limit=500" % (url, self.offset)
        req = Request(url=url, data=data, headers=headers)

        if method:
            req.get_method = lambda: method

        error = False
        try:
            try:
                res = opener.open(req)
            except HTTPError as excpt:
                # It seems that some version of urllib2 are buggy
                if excpt.code <= 299:
                    res = excpt
                else:
                    raise
        except HTTPError as excpt:
            print(repr(excpt))
            error = True
            res = excpt
            # 502: Proxy Error
            # 503: Service Temporarily Unavailable
            if excpt.code >= 500:
                raise RetryError(excpt.code, excpt.read())

        self.retry = int(res.headers.get("Retry-After", self.retry))
        code = res.code
        if code in [201, 202]:
            self.location = res.headers.get("Location", self.location)

        if self.verbose:
            print("Code " + code)
            print("Content-Type " + res.headers.get("Content-Type"))
            print("Content-Length " + res.headers.get("Content-Length"))
            print("Location " + res.headers.get("Location"))

        body = res.read()
        res.close()

        if sys.version_info[0] > 2:
            body = body.decode('utf-8')

        if code in [204]:
            self.last = None
            return None
        else:
            try:
                self.last = json.loads(body)
            except Exception as excpt:
                self.last = {"error": "%s: %s" % (excpt, body)}
                error = True

        if self.verbose:
            print(repr(json.dumps(self.last, indent=4)))

        self.status = self.last.get("status", self.status)

        if self.verbose:
            print("Status " + self.status)

        if "messages" in self.last:
            for n in self.last["messages"]:
                if not self.quiet:
                    print(n)
                self.offset += 1

        if code == 200 and self.status == "complete":
            self.value = self.last
            self.done = True
            if isinstance(self.value, dict) and "result" in self.value:
                self.value = self.value["result"]

        if code in [303]:
            self.value = self.last
            self.done = True

        if "error" in self.last:
            raise APIException("ecmwf.API error 1: %s" % (self.last["error"],))

        if error:
            raise APIException("ecmwf.API error 2: %s" % (res, ))

        return self.last

    def submit(self, url, payload):
        self.call(url, payload, "POST")

    def POST(self, url, payload):
        return self.call(url, payload, "POST")

    def GET(self, url):
        return self.call(url, None, "GET")

    def wait(self):
        if self.verbose:
            print("Sleeping %s second(s)" % (str(self.retry)))
        time.sleep(self.retry)
        self.call(self.location, None, "GET")

    def ready(self):
        return self.done

    def result(self):
        return self.value

    def cleanup(self):
        try:
            if self.location:
                self.call(self.location, None, "DELETE")
        except Exception as e:
            pass


class Ecmwf(object):
    """
    Ecmwf class
    """

    def __init__(self, dates, request):

        # ~> Initialisation
        self.moddates = dates
        self.request = request
        # ~~> inheritence
        self.slf2d = Selafin('')
        self.slf2d.title = ''
        self.slf2d.fole = {}
        self.connection = None
        self.ecmwfdata = None
        self.nx1d = None
        self.ny1d = None
        self.maskx = None
        self.masky = None
        self.nb_freq = None
        self.nb_direct = None
        self.freq = None
        self.dirc = None
        self.typ = None

    def connect_to_ecmwf(self, dataset):

        status = ''
        # ~> Establish connection
        self.connection = Connection(CONFIG['email'], CONFIG['key'],
                                     quiet=True, verbose=False)
        # ~> Verify connection
        user = self.connection.call("%s/%s" % (CONFIG['url'], "who-am-i"))
        print('   ~> access through username: %s\n' %
              (user["full_name"] or "user '%s'" % user["uid"],))
        # ~> Request dataset
        self.connection.submit("%s/%s/requests" % (CONFIG['url'], dataset),
                               self.request)
        status = self.connection.status
        print('   ~> request has been ' + status)
        # ~> Wait for remote processing
        while not self.connection.ready():
            if status != self.connection.status:
                status = self.connection.status
                print('   ~> request remains ' + status + '...')
            self.connection.wait()
        # ~> Request completed
        print('   ~> request is now ' + self.connection.status)
        self.connection.cleanup()

    def download_ecmwf(self):

        result = self.connection.result()
        file_name = self.request.get("target")

        # ~> tries connecting 3 times before stopping
        tries = 0
        while True:

            # ~> downloading file by blocks
            http = urlopen(result["href"])
            f = open(file_name, "wb")
            ibar = 0
            pbar = ProgressBar(maxval=result["size"]).start()
            while True:
                chunk = http.read(1024*1024)
                if not chunk:
                    break
                f.write(chunk)
                ibar += len(chunk)
                pbar.update(ibar)
            f.flush()
            f.close()
            pbar.finish()
            # ~> have I got everything ?
            if ibar == result["size"]:
                break
            if tries == 3:
                raise TelemacException(
                        "    ... exhausted the number "
                        "of download trials.\nYou may wish "
                        "to attempt this again later.")
            print("    ... trying to download the data once more ...")
            tries += 1

    def open_ecmwf(self):

        self.ecmwfdata = netcdf.netcdf_file(self.request.get("target"), 'r')

        # ~~~~ Time records ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ats = self.ecmwfdata.variables['time'][:]-70*365.25*24
        at0 = datetime.fromtimestamp(ats[0]*3600.)
        self.slf2d.datetime = [d for d in at0.timetuple()[0:6]]
        # time record in hours
        self.slf2d.tags = {'times': 3600 * (ats - ats[0])}

    def close_ecmwf(self):
        self.slf2d.fole['hook'].close()

    def set_geometry(self):

        # ~~> 2D grid
        print('   +> set the mesh and connectivity')
        x = self.ecmwfdata.variables['longitude'][:]
        self.nx1d = len(x)
        y = self.ecmwfdata.variables['latitude'][:]
        self.ny1d = len(y)
        # ~~> TODO: a more complicated check should be done here
        #           with a mix of entries
        y_1, x_1, y_2, x_2 = self.request['area'].split('/')

        if float(x_1) < 0.:
            x_1 = float(x_1) + 360.0
        if float(x_2) < 0.:
            x_2 = float(x_2) + 360.0
        self.maskx = np.logical_and(float(x_1) <= x[0], x[0] <= float(x_2))
        if not np.any(self.maskx):
            raise TelemacException(
                    '... your spatial range seems out of bound:\n       '
                    'you asked for [ {}-{}], while x is:\n       {}'
                    '\n\n'.format(x_1, x_2, repr(x)))

        self.masky = np.logical_and(float(y_1) <= y.T[0], y.T[0] <= float(y_2))
        if not np.any(self.masky):
            raise TelemacException(
                    '... your spatial range seems out of bound:\n       '
                    'you asked for [ {}-{}], while x is:\n       {}'
                    '\n\n'.format(y_1, y_2, repr(y)))

        self.slf2d.meshx = np.tile(x, self.ny1d)\
            .reshape(self.ny1d, self.nx1d).T.ravel()
        self.slf2d.meshy = np.tile(y, self.nx1d)

        self.slf2d.nplan = 1
        self.slf2d.ndp2 = 3
        self.slf2d.ndp3 = self.slf2d.ndp2
        self.slf2d.npoin2 = self.nx1d * self.ny1d
        self.slf2d.npoin3 = self.slf2d.npoin2
        self.slf2d.nelem2 = 2*(self.nx1d-1)*(self.ny1d-1)
        self.slf2d.nelem3 = self.slf2d.nelem2
        # ~~> lat,lon correction
        for i in range(self.slf2d.npoin2):
            if self.slf2d.meshx[i] > 180:
                self.slf2d.meshx[i] = self.slf2d.meshx[i] - 360.0
        # for i in range(2172,self.ny1d):
        #   self.slf2d.meshy[i] = 47.0 + ( i-2172 )/18.0

        # ~~> Connectivity
        ielem = 0
        pbar = ProgressBar(maxval=self.slf2d.nelem3).start()
        self.slf2d.ikle3 = np.zeros((self.slf2d.nelem3, self.slf2d.ndp3),
                                    dtype=np.int)
        for i in range(1, self.nx1d):
            for j in range(1, self.ny1d):
                ipoin = (i-1)*self.ny1d + j - 1
                # ~~> first triangle
                self.slf2d.ikle3[ielem][0] = ipoin
                self.slf2d.ikle3[ielem][1] = ipoin + self.ny1d
                self.slf2d.ikle3[ielem][2] = ipoin + 1
                ielem = ielem + 1
                pbar.update(ielem)
                # ~~> second triangle
                self.slf2d.ikle3[ielem][0] = ipoin + self.ny1d
                self.slf2d.ikle3[ielem][1] = ipoin + self.ny1d + 1
                self.slf2d.ikle3[ielem][2] = ipoin + 1
                ielem = ielem + 1
                pbar.update(ielem)
        pbar.finish()

        # ~~> Boundaries
        pbar = ProgressBar(maxval=self.nx1d+self.ny1d).start()
        self.slf2d.ipob3 = np.zeros(self.slf2d.npoin3, dtype=np.int)
        # ~~> along the x-axis (lon)
        for i in range(self.nx1d):
            ipoin = i*self.ny1d
            self.slf2d.ipob3[ipoin] = i + 1
            ipoin = i*self.ny1d - 1
            self.slf2d.ipob3[ipoin] = 2*self.nx1d+(self.ny1d-2) - i
            pbar.update(i)
        # ~~> along the y-axis (alt)
        for i in range(1, self.ny1d):
            ipoin = i
            self.slf2d.ipob3[ipoin] = 2*self.nx1d + 2*(self.ny1d-2) - i + 1
            ipoin = self.ny1d*(self.nx1d-1) + i
            self.slf2d.ipob3[ipoin] = self.nx1d + i
            pbar.update(i+self.nx1d)
        pbar.finish()

        # ~~> Boundary points
        self.slf2d.iparam = [0, 0, 0, 0, 0, 0, 0,
                             2*self.nx1d+2*(self.ny1d-2), 0, 1]

    def put_geometry(self, file_name):

        print('   +> writing up the geometry file')

        self.slf2d.fole = {}
        self.slf2d.fole.update({'hook': open(file_name, 'wb')})
        self.slf2d.fole.update({'name': file_name})
        self.slf2d.fole.update({'endian': ">"})     # big endian
        self.slf2d.fole.update({'float': ('f', 4)})  # single precision

        self.slf2d.varnames = ['RANGE          ']
        self.slf2d.varunits = ['UI             ']
        self.slf2d.nbv1 = len(self.slf2d.varnames)
        self.slf2d.nvar = self.slf2d.nbv1
        self.slf2d.varindex = range(self.slf2d.nvar)

        print('       - Write Selafin geometry')
        self.slf2d.append_header_slf()

        print('       - Write Selafin core')
        varof = self.ecmwfdata.variables['d2fd'].add_offset
        varsf = self.ecmwfdata.variables['d2fd'].scale_factor
        var2d = np.zeros((self.nx1d, self.ny1d), dtype=np.float)
        ibar = 0
        pbar = ProgressBar(maxval=len(self.slf2d.tags['times'])).start()
        for itime in range(len(self.slf2d.tags['times'])):
            self.slf2d.append_core_time_slf(itime)
            var = self.ecmwfdata.variables['d2fd'][itime]
            z = 10 ** (varsf * np.swapaxes(np.swapaxes(var, 1, 3), 0, 2) +
                       varof)
            for j in range(self.ny1d):
                for i in range(self.nx1d):
                    var2d[i, j] = max(z[j][i].ravel())
            self.slf2d.append_core_vars_slf([var2d.ravel()])
            ibar += 1
            pbar.update(ibar)
        pbar.finish()

        self.slf2d.fole['hook'].close()

    def set_spectral(self):

        print('   +> reseting the header of the spectral file')

        print('      - read the spectra definition')
        self.nb_freq = len(self.ecmwfdata.variables['frequency'][:])
        self.nb_direct = len(self.ecmwfdata.variables['direction'][:])
        self.freq = [0.035, 0.038, 0.042, 0.046, 0.051, 0.056, 0.061, 0.067,
                     0.074, 0.081, 0.09, 0.098, 0.108, 0.119, 0.131, 0.144,
                     0.159, 0.174, 0.192, 0.211, 0.232, 0.255, 0.281, 0.309,
                     0.34, 0.374, 0.411, 0.453, 0.498, 0.548]
        #  /!? only for TOMAWC to work
        self.dirc = 7.5 + 15.*np.arange(self.nb_direct) - 7.5

        # ~~> sizes (spectral numbers)
        self.slf2d.nplan = 1
        self.slf2d.ndp2 = 4
        self.slf2d.ndp3 = self.slf2d.ndp2
        self.slf2d.npoin2 = self.nb_direct * self.nb_freq
        self.slf2d.npoin3 = self.slf2d.npoin2
        self.slf2d.nelem2 = self.nb_direct * (self.nb_freq-1)
        self.slf2d.nelem3 = self.slf2d.nelem2
        self.slf2d.nptfr = 2*self.nb_direct
        self.slf2d.iparam = [0, 0, 0, 0, 0, 0, 0, 2*self.nb_direct, 0, 1]

        # ~~> 2D grid (spectral grid) - TODO: use numpy here !
        self.slf2d.meshx = np.zeros(self.slf2d.npoin2, dtype=np.float)
        self.slf2d.meshy = np.zeros(self.slf2d.npoin2, dtype=np.float)
        print('      - set the mesh')
        ipoin = 0
        pbar = ProgressBar(maxval=self.slf2d.npoin2).start()
        for j_f in range(self.nb_freq):
            for j_d in range(self.nb_direct):
                self.slf2d.meshx[j_d+self.nb_direct*j_f] = \
                          self.freq[j_f]*math.sin(math.pi*self.dirc[j_d]/180.)
                self.slf2d.meshy[j_d+self.nb_direct*j_f] = \
                    self.freq[j_f]*math.cos(math.pi*self.dirc[j_d]/180.)
                ipoin += 1
                pbar.update(ipoin)
        pbar.finish()

        # ~~> Connectivity - TODO: use numpy here !
        print('      - set the connectivity')
        ielem = 0
        pbar = ProgressBar(maxval=self.slf2d.nelem3).start()
        self.slf2d.ikle3 = np.zeros((self.slf2d.nelem3, self.slf2d.ndp3),
                                    dtype=np.int)
        for j_f in range(self.nb_freq-1):
            for j_d in range(self.nb_direct):
                self.slf2d.ikle3[ielem][0] = (j_d+1) % self.nb_direct + \
                                                      j_f*self.nb_direct
                ielem += 1
        for ielem in range(self.slf2d.nelem3):
            self.slf2d.ikle3[ielem][1] = ielem
            self.slf2d.ikle3[ielem][2] = ielem + self.nb_direct
            self.slf2d.ikle3[ielem][3] = self.slf2d.ikle3[ielem][0] + \
                self.nb_direct
            pbar.update(ielem)
        pbar.finish()

        # ~~> Boundaries - TODO: use numpy here !
        self.slf2d.ipob3 = np.zeros(self.slf2d.npoin3, dtype=np.int)
        # ~~> along the ?-axis
        for j_d in range(self.nb_direct):
            self.slf2d.ipob3[j_d] = j_d
        for j_d in range(self.nb_direct, 2*self.nb_direct):
            self.slf2d.ipob3[j_d] = self.nb_direct * (self.nb_freq + 1) - j_d

    def append_header_ecmwf(self):

        self.slf2d.varnames = []
        self.slf2d.varunits = []
        # ~~> variables
        self.slf2d.title = ''
        if self.typ == 'wave':
            self.slf2d.varnames = ['WAVE HEIGHT     ',
                                   'WAVE PERIOD     ', 'WAVE DIRECTION  ']
            self.slf2d.varunits = ['M               ',
                                   'S               ', 'DEGREES         ']
        elif self.typ == 'spec':
            for i in range(self.nx1d * self.ny1d):
                self.slf2d.varnames.append(
                          ('F PT '+str(i+1)+'                ')[:16])
                self.slf2d.varunits.append('UI              ')
            print('    - from ', self.slf2d.varnames[0], ' to ',
                  self.slf2d.varnames[-1])
        else:
            self.slf2d.varnames = ['SURFACE PRESSURE',
                                   'WIND VELOCITY U ', 'WIND VELOCITY V ',
                                   'AIR TEMPERATURE ']
            self.slf2d.varunits = ['UI              ',
                                   'M/S             ', 'M/S             ',
                                   'DEGREES         ']
        self.slf2d.nbv1 = len(self.slf2d.varnames)
        self.slf2d.nvar = self.slf2d.nbv1
        self.slf2d.varindex = range(self.slf2d.nvar)

        self.slf2d.append_header_slf()

    def append_core_time_ecmwf(self, itime):
        self.slf2d.append_core_time_slf(itime)

    def append_core_vars_ecmwf(self, itime):
        # Note: this is how you get to the attributes ...
        # ecmwfdata.variables['sp'].ncattrs()
        # in particular ...
        # ecmwfdata.variables['sp'].units
        # ecmwfdata.variables['sp'].missing_value

        if self.typ == 'wave':
            # ~~> WAVE HEIGHT == 'swh'
            var2d = np.swapaxes(self.ecmwfdata.variables['swh'][itime][:],
                                0, 1).ravel()
            varof = self.ecmwfdata.variables['swh'].add_offset
            varsf = self.ecmwfdata.variables['swh'].scale_factor
            self.slf2d.append_core_vars_slf([varsf*var2d+varof])

            # ~~> SIGNIFICANT WAVE PERIOD == 'mwp'
            var2d = np.swapaxes(self.ecmwfdata.variables['mwp'][itime][:],
                                0, 1).ravel()
            varof = self.ecmwfdata.variables['mwp'].add_offset
            varsf = self.ecmwfdata.variables['mwp'].scale_factor
            self.slf2d.append_core_vars_slf([varsf*var2d+varof])

            # ~~> MEAN WAVE DIRECTION == 'mwd'
            var2d = np.swapaxes(self.ecmwfdata.variables['mwd'][itime][:],
                                0, 1).ravel()
            varof = self.ecmwfdata.variables['mwd'].add_offset
            varsf = self.ecmwfdata.variables['mwd'].scale_factor
            self.slf2d.append_core_vars_slf([varsf*var2d+varof])

        elif self.typ == 'spec':
            varof = self.ecmwfdata.variables['d2fd'].add_offset
            varsf = self.ecmwfdata.variables['d2fd'].scale_factor
            var = self.ecmwfdata.variables['d2fd'][itime]
            z = 10 ** (varsf * np.swapaxes(np.swapaxes(var, 1, 3), 0, 2) +
                       varof)
            for j in range(self.ny1d):
                for i in range(self.nx1d):
                    self.slf2d.append_core_vars_slf([z[j][i].ravel()])

        else:
            # ~~> SURFACE PRESSURE == 'sp'
            var2d = np.swapaxes(self.ecmwfdata.variables['sp'][itime][:],
                                0, 1).ravel()
            varof = self.ecmwfdata.variables['sp'].add_offset
            varsf = self.ecmwfdata.variables['sp'].scale_factor
            # print( self.ecmwfdata.variables['sp'].units )
            self.slf2d.append_core_vars_slf([varsf*var2d+varof])

            # ~~> WIND VELOCITY U == 'u10'
            var2d = np.swapaxes(self.ecmwfdata.variables['u10'][itime][:],
                                0, 1).ravel()
            varof = self.ecmwfdata.variables['u10'].add_offset
            varsf = self.ecmwfdata.variables['u10'].scale_factor
            # print( self.ecmwfdata.variables['u10'].units )
            self.slf2d.append_core_vars_slf([varsf*var2d+varof])

            # ~~> WIND VELOCITY V == 'v10'
            var2d = np.swapaxes(self.ecmwfdata.variables['v10'][itime][:],
                                0, 1).ravel()
            varof = self.ecmwfdata.variables['v10'].add_offset
            varsf = self.ecmwfdata.variables['v10'].scale_factor
            # print( self.ecmwfdata.variables['v10'].units )
            self.slf2d.append_core_vars_slf([varsf*var2d+varof])

            # ~~> AIR TEMPERATURE == 't2m'
            var2d = np.swapaxes(self.ecmwfdata.variables['t2m'][itime][:],
                                0, 1).ravel()
            varof = self.ecmwfdata.variables['t2m'].add_offset
            varsf = self.ecmwfdata.variables['t2m'].scale_factor
            # Kelvin to Celsius
            self.slf2d.append_core_vars_slf([varsf*var2d+varof-273.15])

    def put_content(self, file_name, stream, showbar=True):

        # ~~> netcdf reader
        self.typ = stream

        if self.typ == 'spec':
            # ~~> new Selafin writer
            self.slf2d.fole = {}
            self.slf2d.fole.update({'hook': open(file_name, 'wb')})
            self.slf2d.fole.update({'name': file_name})
            self.slf2d.fole.update({'endian': ">"})     # big endian
            self.slf2d.fole.update({'float': ('f', 4)})  # single precision

            print('     +> Write Selafin header')
            self.append_header_ecmwf()

            print('     +> Write Selafin core')
            ibar = 0
            if showbar:
                pbar = ProgressBar(maxval=len(self.slf2d.tags['times']))\
                    .start()
            for itime in range(len(self.slf2d.tags['times'])):
                self.append_core_time_ecmwf(itime)
                self.append_core_vars_ecmwf(ibar)
                ibar += 1
                if showbar:
                    pbar.update(ibar)
            self.slf2d.fole['hook'].close()
            if showbar:
                pbar.finish()

        else:
            # ~~> new Selafin writer
            self.slf2d.fole = {}
            self.slf2d.fole.update({'hook': open(file_name, 'wb')})
            self.slf2d.fole.update({'name': file_name})
            self.slf2d.fole.update({'endian': ">"})     # big endian
            self.slf2d.fole.update({'float': ('f', 4)})  # single precision

            print('     +> Write Selafin header')
            self.append_header_ecmwf()

            print('     +> Write Selafin core')
            ibar = 0
            if showbar:
                pbar = ProgressBar(maxval=len(self.slf2d.tags['times']))\
                    .start()
            for itime in range(len(self.slf2d.tags['times'])):
                self.append_core_time_ecmwf(itime)
                self.append_core_vars_ecmwf(ibar)
                ibar += 1
                if showbar:
                    pbar.update(ibar)
            self.slf2d.fole['hook'].close()
            if showbar:
                pbar.finish()
