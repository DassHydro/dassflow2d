"""@author Tobias Bieniek
        Tobias.Bieniek@gmx.de
        https://github.com/Turbo87/utm

    @brief
        Bidirectional UTM-WGS84 converter for python

    @details
    > Ellipsoid name,    Equatorial Radius,  square of eccentricity
    "Airy",                    R =  6377563, E = 0.00667054,
    "Australian National",     R =  6378160, E = 0.006694542,
    "Bessel 1841",             R =  6377397, E = 0.006674372,
    "Bessel 1841 (Nambia] ",   R =  6377484, E = 0.006674372,
    "Clarke 1866",             R =  6378206, E = 0.006768658,
    "Clarke 1880",             R =  6378249, E = 0.006803511,
    "Everest",                 R =  6377276, E = 0.006637847,
    "Fischer 1960 (Mercury] ", R =  6378166, E = 0.006693422,
    "Fischer 1968",            R =  6378150, E = 0.006693422,
    "GRS 1967",                R =  6378160, E = 0.006694605,
    "GRS 1980",                R =  6378137, E = 0.00669438,
    "Helmert 1906",            R =  6378200, E = 0.006693422,
    "Hough",                   R =  6378270, E = 0.00672267,
    "International",           R =  6378388, E = 0.00672267,
    "Krassovsky",              R =  6378245, E = 0.006693422,
    "Modified Airy",           R =  6377340, E = 0.00667054,
    "Modified Everest",        R =  6377304, E = 0.006637847,
    "Modified Fischer 1960",   R =  6378155, E = 0.006693422,
    "South American 1969",     R =  6378160, E = 0.006694542,
    "WGS 60",                  R =  6378165, E = 0.006693422,
    "WGS 66",                  R =  6378145, E = 0.006694542,
    "WGS-72",                  R =  6378135, E = 0.006694318,
    "WGS-84",                  R =  6378137, E = 0.00669438
"""
# Copyright (C) 2012 Tobias Bieniek <Tobias.Bieniek@gmx.de>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#    The above copyright notice and this permission notice shall be included in
#    all copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#    DEALINGS IN THE SOFTWARE.
#
# For most use cases in this module, numpy is indistinguishable
# from math, except it also works on numpy arrays
import numpy as np
from utils.exceptions import TelemacException

__all__ = ['to_latlon', 'from_latlon']

K0 = 0.9996

E = 0.00669438
E2 = E * E
E3 = E2 * E
E_P2 = E / (1.0 - E)

SQRT_E = np.sqrt(1 - E)
_E = (1 - SQRT_E) / (1 + SQRT_E)
_E2 = _E * _E
_E3 = _E2 * _E
_E4 = _E3 * _E
_E5 = _E4 * _E

M1 = (1 - E / 4 - 3 * E2 / 64 - 5 * E3 / 256)
M2 = (3 * E / 8 + 3 * E2 / 32 + 45 * E3 / 1024)
M3 = (15 * E2 / 256 + 45 * E3 / 1024)
M4 = (35 * E3 / 3072)

P2 = (3. / 2 * _E - 27. / 32 * _E3 + 269. / 512 * _E5)
P3 = (21. / 16 * _E2 - 55. / 32 * _E4)
P4 = (151. / 96 * _E3 - 417. / 128 * _E5)
P5 = (1097. / 512 * _E4)

R = 6378137

ZONE_LETTERS = "CDEFGHJKLMNPQRSTUVWXX"


def in_bounds(x, lower, upper, upper_strict=False):
    """
    Check that x is in bound of [lower, upper]

    @param x (float/int) value to check
    @param lower (float/int) Min of bound
    @param upper (float/int) Max of bound
    @param upper_strict (bool) Max included or not

    @return (bool)
    """
    if upper_strict:
        return lower <= np.min(x) and np.max(x) < upper
    return lower <= np.min(x) and np.max(x) <= upper


def check_valid_zone(zone_number, zone_letter):
    """
    Check that zone number and zone letter are valid

    @param zone_number (int) Zone number
    @param zone_letter (str) Zone letter

    @returns (bool) True if valis
    """
    if not 1 <= zone_number <= 60:
        raise TelemacException(
            'zone number out of range (must be between 1 and 60)')

    if zone_letter:
        zone_letter = zone_letter.upper()

        if not 'C' <= zone_letter <= 'X' or zone_letter in ['I', 'O']:
            raise TelemacException(
                'zone letter out of range (must be between C and X)')


def mixed_signs(x):
    """
    Check if x contains mixed signs

    @param x (np.array) The array to check

    @returns (bool) True id mixed signs
    """
    return np.min(x) < 0 and np.max(x) >= 0


def negative(x):
    """
    Returns True if x contains only negative values

    @param x (np.array) x

    @returns (bool)
    """
    return np.max(x) < 0


def to_latlon(easting, northing, zone_number, zone_letter=None,
              northern=None, strict=True):
    """
    This function convert an UTM coordinate into Latitude and Longitude

    @param easting (int)
            Easting value of UTM coordinate

    @param northing (int)
            Northing value of UTM coordinate

    @param zone number (int)
            Zone Number is represented with global map numbers of an UTM Zone
            Numbers Map. More information see utmzones [1]_

    @param zone_letter: str
            Zone Letter can be represented as string values. Where UTM Zone
            Designators can be accessed in [1]_

    @param northern: bool
            You can set True or False to set this parameter. Default is None


       .. _[1]: http://www.jaworski.ca/utmzones.htm

    @returns (np.array, np.array) longitude, latitude
    """
    if not zone_letter and northern is None:
        raise ValueError('either zone_letter or northern needs to be set')

    if zone_letter and northern is not None:
        raise ValueError('set either zone_letter or northern, but not both')

    if strict:
        if not in_bounds(easting, 100000, 1000000, upper_strict=True):
            raise TelemacException(
                'easting out of range (must be between 100.000m and 999.999m)')
        if not in_bounds(northing, 0, 10000000):
            raise TelemacException(
                'northing out of range (must be between 0m and 10.000.000m)')

    check_valid_zone(zone_number, zone_letter)

    if zone_letter:
        zone_letter = zone_letter.upper()
        northern = (zone_letter >= 'N')

    x = easting - 500000
    y = northing

    if not northern:
        y -= 10000000

    m = y / K0
    m_u = m / (R * M1)

    p_rad = (m_u +
             P2 * np.sin(2 * m_u) +
             P3 * np.sin(4 * m_u) +
             P4 * np.sin(6 * m_u) +
             P5 * np.sin(8 * m_u))

    p_sin = np.sin(p_rad)
    p_sin2 = p_sin * p_sin

    p_cos = np.cos(p_rad)

    p_tan = p_sin / p_cos
    p_tan2 = p_tan * p_tan
    p_tan4 = p_tan2 * p_tan2

    ep_sin = 1 - E * p_sin2
    ep_sin_sqrt = np.sqrt(1 - E * p_sin2)

    n = R / ep_sin_sqrt
    r_0 = (1 - E) / ep_sin

    c_0 = _E * p_cos**2
    c_2 = c_0 * c_0

    d_0 = x / (n * K0)
    d_2 = d_0 * d_0
    d_3 = d_2 * d_0
    d_4 = d_3 * d_0
    d_5 = d_4 * d_0
    d_6 = d_5 * d_0

    latitude = (p_rad - (p_tan / r_0) *
                (d_2 / 2 -
                 d_4 / 24 * (5 + 3 * p_tan2 + 10 * c_0 - 4 * c_2 - 9 * E_P2))
                + d_6 / 720 * (61 + 90 * p_tan2 + 298 * c_0 + 45 * p_tan4
                               - 252 * E_P2 - 3 * c_2))

    longitude = (d_0 -
                 d_3 / 6 * (1 + 2 * p_tan2 + c_0) +
                 d_5 / 120 * (5 - 2 * c_0 + 28 * p_tan2 - 3 * c_2 + 8 * E_P2
                              + 24 * p_tan4)) / p_cos

    return (np.degrees(longitude)
            + zone_number_to_central_longitude(zone_number),
            np.degrees(latitude))


def from_latlon(latitude, longitude, force_zone_number=None,
                force_zone_letter=None):
    """
    This function convert Latitude and Longitude to UTM coordinate

    @param latitude (np.array)
             Latitude between 80 deg S and 84 deg N, e.g. (-80.0 to 84.0)

    @param longitude (np.array)
             Longitude between 180 deg W and 180 deg E, e.g. (-180.0 to 180.0).

    @param force_zone_number (int)
             Zone Number is represented with global map numbers of an UTM Zone
             Numbers Map. You may force conversion including one UTM Zone
             Number.
             More information see utmzones [1]_
    @param force_zone_letter (str)
             Zone letter to be forced

       .. _[1]: http://www.jaworski.ca/utmzones.htm

    @returns (np.array, np.array, int) easting, northing, zone
    """
    if not in_bounds(latitude, -80.0, 84.0):
        raise TelemacException(
            'latitude out of range (must be between 80 deg S and 84 deg N)')
    if not in_bounds(longitude, -180.0, 180.0):
        raise TelemacException(
            'longitude out of range (must be between 180 deg W and 180 deg E)')
    if force_zone_number is not None:
        check_valid_zone(force_zone_number, force_zone_letter)

    lat_rad = np.radians(latitude)
    lat_sin = np.sin(lat_rad)
    lat_cos = np.cos(lat_rad)

    lat_tan = lat_sin / lat_cos
    lat_tan2 = lat_tan * lat_tan
    lat_tan4 = lat_tan2 * lat_tan2

    if force_zone_number is None:
        zone_number = latlon_to_zone_number(latitude, longitude)
    else:
        zone_number = force_zone_number

    if force_zone_letter is None:
        zone_letter = latitude_to_zone_letter(latitude)
    else:
        zone_letter = force_zone_letter

    lon_rad = np.radians(longitude)
    central_lon = zone_number_to_central_longitude(zone_number)
    central_lon_rad = np.radians(central_lon)

    n = R / np.sqrt(1 - E * lat_sin**2)
    c_0 = E_P2 * lat_cos**2

    a_0 = lat_cos * (lon_rad - central_lon_rad)
    a_2 = a_0 * a_0
    a_3 = a_2 * a_0
    a_4 = a_3 * a_0
    a_5 = a_4 * a_0
    a_6 = a_5 * a_0

    m = R * (M1 * lat_rad -
             M2 * np.sin(2 * lat_rad) +
             M3 * np.sin(4 * lat_rad) -
             M4 * np.sin(6 * lat_rad))

    easting = K0 * n * (a_0 +
                        a_3 / 6 * (1 - lat_tan2 + c_0) +
                        a_5 / 120 * (5 - 18*lat_tan2 + lat_tan4
                                     + 72*c_0 - 58*E_P2)) + 500000

    northing = K0 * (m + n*lat_tan*(a_2 / 2 +
                                    a_4 / 24 * (5 - lat_tan2 + 9*c_0
                                                + 4*c_0**2) +
                                    a_6 / 720 * (61 - 58*lat_tan2 + lat_tan4
                                                 + 600*c_0 - 330 * E_P2)))

    if mixed_signs(latitude):
        raise ValueError("latitudes must all have the same sign")

    if negative(latitude):
        northing += 10000000

    return easting, northing, zone_number, zone_letter


def latitude_to_zone_letter(latitude):
    """
    Identify zone letter from lattitude

    @param latitude (np.array) latitude

    @returns (str) zone letter
    """
    # If the input is a numpy array, just use the first element
    # User responsibility to make sure that all points are in one zone
    latitude = latitude.flat[0]

    if -80 <= latitude <= 84:
        return ZONE_LETTERS[int(latitude + 80) >> 3]

    return None


def latlon_to_zone_number(latitude, longitude):
    """
    Identify zone number from lattitude and longitude

    @param latitude (np.array) latitude
    @param longitude (np.array) longitude

    @returns (int) zone number
    """
    # If the input is a numpy array, just use the first element
    # User responsibility to make sure that all points are in one zone
    latitude = latitude.flat[0]
    longitude = longitude.flat[0]

    if 56 <= latitude < 64 and 3 <= longitude < 12:
        return 32

    if 72 <= latitude <= 84 and longitude >= 0:
        if longitude < 9:
            return 31
        if longitude < 21:
            return 33
        if longitude < 33:
            return 35
        if longitude < 42:
            return 37

    return int((longitude + 180) / 6) + 1


def zone_number_to_central_longitude(zone_number):
    """
    Returns central longitude of a zone number

    @param zone_number (int) zone number

    @returns (int) central longitude
    """
    return (zone_number - 1) * 6 - 180 + 3
