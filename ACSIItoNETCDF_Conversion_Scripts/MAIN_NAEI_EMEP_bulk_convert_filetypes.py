#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:46:59 2019

Script for creating netcdf files from the NAEI text files, 
converting from BNG to lat / lon grids in the process.


@author: mbessdl2
"""
#import sys
import pandas as pd
import numpy as np
import xarray as xr
#from bng_to_latlon import OSGB36toWGS84

from math import sqrt, pi, sin, cos, tan, atan2
from numba import jit


#%% general data settings
root_dir = "../../NAEI_2016/"
root_tail = "16.asc"

species_dir = ["CH4","CO","HCl","NH3","NMVOC","NOx","PM01","PM1","PM25","PM10","SO2"]
species_str = {"CH4"  :"ch4",\
               "CO"   :"co",\
               "HCl"  :"hcl",\
               "NH3"  :"nh3",\
               "NMVOC":"voc",\
               "NOx"  :"nox",\
               "PM01" :"pm0_1",\
               "PM1"  :"pm1",\
               "PM25" :"pm2_5",\
               "PM10" :"pm10",\
               "SO2"  :"so2"}
netcdf_names ={"CH4"  :"ch4",\
               "CO"   :"co",\
               "HCl"  :"hcl",\
               "NH3"  :"nh3",\
               "NMVOC":"nmvoc",\
               "NOx"  :"nox",\
               "PM01" :"pm0_1",\
               "PM1"  :"pm1",\
               "PM25" :"pm2_5",\
               "PM10" :"pm10",\
               "SO2"  :"so2"}

input_paths = {value: root_dir+value+"/" for value in species_dir}
input_tails = {value: species_str[value]+root_tail for value in species_dir}

output_paths = input_paths



header_end = "NODATA"

var_names = ["agric","domcom","energyprod","indcom","indproc","nature",\
             "offshore","othertrans","roadtrans","solvents","total",\
             "totarea","waste"]



#%% copy of the bng_to_latlon function, with numba optimisation
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



#%% functions 

#
#  Function for extracting the header data from an ASCII file, and
#     also to count how deep the header is.
#
#  This function requires a string value to be passed, which defines
#     the string starts the final line of the header.
#
def read_header(filepath,filename,headend):

    head_dict = {}
    
	# loop through header, splitting by white space, and saving the 
	#   header information as attributes.
	# Note: we assume that each header line has two values, and that the 2nd is numeric
    with open (filepath+filename, 'rt') as myfile:
        for myline in myfile:
            head_data = myline.split()
            if len(head_data) == 2:
                head_dict[head_data[0]] = int(head_data[1])
            else:
                break
            
    return(head_dict)


#
#  Function for extracting table of data from ASCII file, and converting to an
#  xarray data array
#
#def read_table(filepath,filename,datatype,headers):
#
#    data_table = pd.read_csv(filepath+filename,header=None,\
#                        skiprows=len(headers),delim_whitespace=True,na_values=headers['NODATA_value'])




def calculate_lat_lon_grids(headers):

    
    # set up the northing & easting data arrays
    e_data = xr.DataArray(np.zeros((headers['nrows'],headers['ncols'])),dims=('lat','lon'))
    n_data = xr.DataArray(np.zeros((headers['nrows'],headers['ncols'])),dims=('lat','lon'))
    
    latitude = xr.DataArray(np.zeros((headers['nrows'],headers['ncols'])),dims=('lat','lon'))
    longitude = xr.DataArray(np.zeros((headers['nrows'],headers['ncols'])),dims=('lat','lon'))
    
    
    # create 1D easting and westing arrays, switching from lower-left position to the middle of gridcell
    eastings  = xr.DataArray( \
                    np.arange(headers['xllcorner'],\
                          (headers['xllcorner']+headers['ncols']*headers['cellsize']),\
                          headers['cellsize']) \
                              + int(headers['cellsize']/2) \
                    ,dims=('lon'))
    northings = xr.DataArray( \
                    np.arange(headers['yllcorner'],\
                          (headers['yllcorner']+headers['nrows']*headers['cellsize']),\
                          headers['cellsize']) \
                              + int(headers['cellsize']/2) \
                    ,dims=('lat'))
    
   
    
    # copy the 1D arrays into the 2D arrays
    for ii in np.arange(0,headers['nrows']):
        e_data[ii,:] = eastings
    for ii in np.arange(0,headers['ncols']):
        n_data[:,ii] = northings

    
    for ii in np.arange(0,headers['nrows']):
        for jj in np.arange(0,headers['ncols']):
            (lat2,lon2) = OSGB36toWGS84(e_data.values[ii,jj],n_data.values[ii,jj])            
            #print(ii,jj,lat2,lon2)
            latitude[ii,jj]  = lat2
            longitude[ii,jj] = lon2
    

    return(latitude,longitude,e_data,n_data)





#%% read the header data, and create the latitude, longitude, eastings and northings grids

## load the header information (from the first file, and assume
##              all data files for all species have the same headers!) 
headers = read_header(input_paths['CH4'],var_names[0]+input_tails['CH4'],header_end)

(lat4,lon4,east4,north4) = calculate_lat_lon_grids(headers)

#%% start creation of the dataset which will form the final files

# template for emission data arrays (with time dimension)
emiss_template = xr.DataArray(np.zeros((2,headers['nrows'],headers['ncols'])),dims=('time','lat','lon'))

# date template array
date = xr.DataArray([19000101, 21000101],dims=('time'))


# create template dataset
ds_template = xr.Dataset({'latitude':lat4,'longitude':lon4,\
                          'eastings':east4,'northings':north4,\
                          'date':date})


#%% processing the data files, and outputting netcdf files

for species in species_dir:
    
    # copy dataset template
    ds_working = ds_template.copy(deep=True)

    
    # process each emission dataset for this species
    for emiss in var_names:

        test_file = input_paths[species]+emiss+input_tails[species]
        emiss_data = emiss_template.copy(deep=True)

        try:
            test_data = pd.read_csv(test_file,header=None,\
                        skiprows=len(headers),delim_whitespace=True,na_values=headers['NODATA_value'])
            print('opened '+emiss+' for '+species)
            emiss_data[0,:,:] = test_data.values[::-1,:] # flip the emission data vertically, 
            emiss_data[1,:,:] = test_data.values[::-1,:] #          to match the lat/lon grid
        except:
            print('skipping '+emiss+' for '+species)

        ds_working[emiss] = emiss_data

    # write to netcdf file
    outfile = output_paths[species]+netcdf_names[species]+'_data.nc'
    ds_working.to_netcdf(path=outfile,mode='w')



