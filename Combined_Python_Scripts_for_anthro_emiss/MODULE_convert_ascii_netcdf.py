#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur May 23 2019

Functions for converting NAEI data from ASCII to NETCDF

@author: mbessdl2
"""

from MODULE_bng_to_latlon_local import OSGB36toWGS84  

import xarray as xr
import numpy as np
import pandas as pd

from sys import exit


def read_header(species_dir,input_paths,var_names,input_tails,header_end):
    """
    Function for extracting the header data from an ASCII file, and
    also to count how deep the header is.

    This function requires a string value to be passed, which defines
    the string starts the final line of the header.
    """

    head_dict = {}
    
    for species in species_dir:
        for emiss in var_names:
            input_file = input_paths[species]+emiss+input_tails[species]
            
            try:    
                # loop through header, splitting by white space, and saving the 
                #   header information as attributes.
                # Note: we assume that each header line has two values, and that the 2nd is numeric
                with open (input_file, 'rt') as myfile:
                    for myline in myfile:
                        head_data = myline.split()
                        if len(head_data) == 2:
                            head_dict[head_data[0]] = int(head_data[1])
                        else:
                            break
                return(head_dict)
            except:
                pass
    
    print('no input files found, failed at reading headers')
    exit(1)


def create_emission_var_template(headers,Atts):
    """
        Function for creating an emission variable template.
    """

    # template for emission data arrays (with time dimension)
    emiss_template = xr.DataArray(np.zeros((2,headers['nrows'],headers['ncols'])),dims=('time','lat','lon'),attrs=Atts)

    return(emiss_template)




def calculate_lat_lon_grids(headers):
    """
    Function for creating the latitude and longitude grids from provided NAEI header data.

    This requires a copy of OSGB36toWGS84 - use the local copy optimised with numba
    
    """
    
    latitude = xr.DataArray(np.zeros((headers['nrows'],headers['ncols'])),dims=('lat','lon'))
    longitude = xr.DataArray(np.zeros((headers['nrows'],headers['ncols'])),dims=('lat','lon'))
    
    latitude.attrs['long_name'] = 'Latitude'
    latitude.attrs['units'] = 'degrees_north'
    longitude.attrs['long_name'] = 'Longitude'
    longitude.attrs['units'] = 'degrees_east'
    
    
    # create 1D easting and westing arrays, switching from lower-left position to the middle of gridcell
    eastings  = xr.DataArray( \
                    np.arange(headers['xllcorner'],\
                          (headers['xllcorner']+headers['ncols']*headers['cellsize']),\
                          headers['cellsize']) \
                              + headers['cellsize']*0.5 \
                    ,dims=('lon'))
    northings = xr.DataArray( \
                    np.arange(headers['yllcorner'],\
                          (headers['yllcorner']+headers['nrows']*headers['cellsize']),\
                          headers['cellsize']) \
                              + headers['cellsize']*0.5 \
                    ,dims=('lat'))
    
    eastings.attrs['long_name'] = 'Eastings'
    eastings.attrs['units'] = 'metres east from origin'
    northings.attrs['long_name'] = 'Northings'
    northings.attrs['units'] = 'metres north from origin'
    
    
    
    # copy the 1D arrays into the 2D arrays
    #for ii in np.arange(0,headers['nrows']):
    #    e_data[ii,:] = eastings
    #for ii in np.arange(0,headers['ncols']):
    #    n_data[:,ii] = northings

    
    for ii in np.arange(0,headers['nrows']):
        print('calculating lat/lon for row ',ii)
        for jj in np.arange(0,headers['ncols']):
            (lat2,lon2) = OSGB36toWGS84(eastings.values[jj],northings.values[ii])            
            #print(ii,jj,lat2,lon2)
            latitude[ii,jj]  = lat2
            longitude[ii,jj] = lon2
    

    return(latitude,longitude,eastings,northings)




def create_NAEI_data_template(headers,Global_Atts):
    """
        Function for creating the NAEI data template.
    
        Uses the calculate_lat_lon_grids function.
    """

    (lat4,lon4,east4,north4) = calculate_lat_lon_grids(headers)

    #%% start creation of the dataset which will form the final files

    # date template array
    date = xr.DataArray([19000101, 21000101],dims=('time'))
    date.attrs['long_name'] = 'Date'
    date.attrs['units'] = 'YYYYMMDD'


    # create template dataset
    ds_template = xr.Dataset({'latitude':lat4,'longitude':lon4,\
                              'eastings':east4,'northings':north4,\
                              'date':date},attrs=Global_Atts)
                              
    return(ds_template)



def load_emission_data(species_dir,species_str,input_paths,var_names,input_tails,\
                                            header_end,Atts,Atts_Sector,Global_Atts,\
                                            emep_names):
    """
    Creates emission and data set templates, loads emission data into 
    these, and returns a dictionary containing all the xarray datasets.
    """

    ## read the header information, first finding a file that exists
    headers = read_header(species_dir,input_paths,var_names,input_tails,header_end)


    ## create the emission variable template (with appropriate attributes)
    emiss_template = create_emission_var_template(headers,Atts)


    ## load the basic file template (containing north/east and lat/lon data)
    data_template = create_NAEI_data_template(headers,Global_Atts)


    ## create dataset dictionary
    all_datasets = {}
    
    
    
    # step through the species in our list, creating a dataset for each one
    for species in species_dir:
    
        # copy dataset template
        all_datasets[species] = data_template.copy(deep=True)

    
        # process each emission dataset for this species
        for emiss in var_names:

            test_file = input_paths[species]+emiss+input_tails[species]
            emiss_data = emiss_template.copy(deep=True)
            
            # add the EMEP specific sector code and species name (where needed)
            try:
                emiss_data.attrs['sector'] = Atts_Sector[emiss]
                emiss_data.attrs['species'] = emep_names[species]
            except:
                pass

            try:
                test_data = pd.read_csv(test_file,header=None,\
                            skiprows=len(headers),delim_whitespace=True,na_values=headers['NODATA_value'])
                test_data = test_data.fillna(0.0)
                print('opened '+emiss+' for '+species)
                emiss_data[0,:,:] = test_data.values[::-1,:] # flip the emission data vertically, 
                emiss_data[1,:,:] = test_data.values[::-1,:] #          to match the lat/lon grid
            except:
                print('skipping '+emiss+' for '+species)

            all_datasets[species][emiss] = emiss_data


    return(all_datasets)




