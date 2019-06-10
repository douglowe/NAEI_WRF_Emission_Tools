#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur May 23 2019

Local copy of bng_to_latlon, to allow for use of numpa functions.

@author: mbessdl2
"""


from math import sqrt, pi, sin, cos, tan, atan2
from numba import jit


#
#  see https://github.com/fmalina/bng_latlon for original version
#

@jit(nopython=True)
def OSGB36toWGS84(E, N):
    """ Accept The Ordnance Survey National Grid eastings and northings.
    Return latitude and longitude coordinates.

    Usage:
    >>> from bng_to_latlon import OSGB36toWGS84
    >>> OSGB36toWGS84(538890, 177320)
    (51.47779538331092, -0.0014016837826672265)
    >>> OSGB36toWGS84(352500.2, 401400)
    (53.507129843104195, -2.7176599627343263)
    """
    # The Airy 1830 semi-major and semi-minor axes used for OSGB36 (m)
    a, b = 6377563.396, 6356256.909
    F0 = 0.9996012717  # scale factor on the central meridian

    # Latitude and longtitude of true origin (radians)
    lat0 = 49.0*pi/180.0
    lon0 = -2.0*pi/180.0  # longtitude of central meridian

    # Northing & easting of true origin (m)
    N0, E0 = -100000.0, 400000.0
    e2 = 1.0 - (b*b)/(a*a)  # eccentricity squared
    n = (a-b)/(a+b)

    # Initialise the iterative variables
    lat, M = lat0, 0.0

    while N-N0-M >= 0.00001:  # Accurate to 0.01mm
        lat = (N-N0-M)/(a*F0) + lat
        M1 = (1.0 + n + (5.0/4.0)*n**2.0 + (5.0/4.0)*n**3.0) * (lat-lat0)
        M2 = (3.0*n + 3.0*n**2.0 + (21.0/8.0)*n**3.0) * sin(lat-lat0) * cos(lat+lat0)
        M3 = ((15.0/8.0)*n**2.0 + (15.0/8.0)*n**3.0) * sin(2.0*(lat-lat0)) * cos(2*(lat+lat0))
        M4 = (35.0/24.0)*n**3.0 * sin(3.0*(lat-lat0)) * cos(3.0*(lat+lat0))
        # meridional arc
        M = b * F0 * (M1 - M2 + M3 - M4)

    # transverse radius of curvature
    nu = a*F0/sqrt(1.0-e2*sin(lat)**2.0)

    # meridional radius of curvature
    rho = a*F0*(1.0-e2)*(1-e2*sin(lat)**2.0)**(-1.5)
    eta2 = nu/rho-1.0

    sec_lat = 1.0/cos(lat)
    VII = tan(lat)/(2.0*rho*nu)
    VIII = tan(lat)/(24.0*rho*nu**3)*(5.0+3.0*tan(lat)**2.0+eta2-9.0*tan(lat)**2.0*eta2)
    IX = tan(lat)/(720.0*rho*nu**5.0)*(61.0+90.0*tan(lat)**2.0+45.0*tan(lat)**4.0)
    X = sec_lat/nu
    XI = sec_lat/(6.0*nu**3.0)*(nu/rho+2.0*tan(lat)**2.0)
    XII = sec_lat/(120.0*nu**5.0)*(5.0+28.0*tan(lat)**2.0+24.0*tan(lat)**4.0)
    XIIA = sec_lat/(5040.0*nu**7.0)*(61.0+662.0*tan(lat)**2.0+1320.0*tan(lat)**4.0+720.0*tan(lat)**6.0)
    dE = E-E0

    # These are on the wrong ellipsoid currently: Airy 1830 (denoted by _1)
    lat_1 = lat - VII*dE**2.0 + VIII*dE**4.0 - IX*dE**6.0
    lon_1 = lon0 + X*dE - XI*dE**3.0 + XII*dE**5.0 - XIIA*dE**7.0

    # Want to convert to the GRS80 ellipsoid.
    # First convert to cartesian from spherical polar coordinates
    H = 0.0  # Third spherical coord.
    x_1 = (nu/F0 + H)*cos(lat_1)*cos(lon_1)
    y_1 = (nu/F0 + H)*cos(lat_1)*sin(lon_1)
    z_1 = ((1.0-e2)*nu/F0 + H)*sin(lat_1)

    # Perform Helmut transform (to go between Airy 1830 (_1) and GRS80 (_2))
    s = -20.4894*10.0**-6.0  # The scale factor -1
    # The translations along x, y, z axes respectively
    tx, ty, tz = 446.448, -125.157, + 542.060
    # The rotations along x, y, z respectively (in seconds)
    rxs, rys, rzs = 0.1502, 0.2470, 0.8421

    # convert seconds to radians
    def sec_to_rad(x): return x*pi/(180.0*3600.0)

    rx, ry, rz = [sec_to_rad(x) for x in (rxs, rys, rzs)]  # (in radians)
    x_2 = tx + (1.0+s)*x_1 + (-rz)*y_1 + (ry)*z_1
    y_2 = ty + (rz)*x_1 + (1.0+s)*y_1 + (-rx)*z_1
    z_2 = tz + (-ry)*x_1 + (rx)*y_1 + (1.0+s)*z_1

    # Back to spherical polar coordinates from cartesian
    # Need some of the characteristics of the new ellipsoid

    # The GSR80 semi-major and semi-minor axes used for WGS84(m)
    a_2, b_2 = 6378137.000, 6356752.3141
    e2_2 = 1.0 - (b_2*b_2)/(a_2*a_2)  # The eccentricity of the GRS80 ellipsoid
    p = sqrt(x_2**2.0 + y_2**2.0)

    # Lat is obtained by an iterative proceedure:
    lat = atan2(z_2, (p*(1.0-e2_2)))  # Initial value
    latold = 2.0*pi
    while abs(lat - latold) > 10.0**-16.0:
        lat, latold = latold, lat
        nu_2 = a_2/sqrt(1.0-e2_2*sin(latold)**2.0)
        lat = atan2(z_2+e2_2*nu_2*sin(latold), p)

    # Lon and height are then pretty easy
    lon = atan2(y_2, x_2)
    H = p/cos(lat) - nu_2

    # Uncomment this line if you want to print the results
    # print([(lat-lat_1)*180/pi, (lon - lon_1)*180/pi])

    # Convert to degrees
    lat = lat*180.0/pi
    lon = lon*180.0/pi

    # Job's a good'n.
    return lat, lon

