#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:46:59 2019

Script for creating EMEP compliant netcdf files from 
the NAEI text files, converting from BNG to lat / lon 
grids in the process, and loading the point source data.


@author: mbessdl2
"""

#import sys
#sys.path.append('./')

#from pathlib import Path
#import os

import MODULE_convert_ascii_netcdf as can
import MODULE_point_source_apportionment as psa

import numpy as np
import xarray as xr
#import pandas as pd





## configuration information

root_dir = "/mnt/iusers01/support/mbessdl2/Projects/Emissions_Processing/Emission_Datasets/NAEI_2016/"
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
               "PM25" :"pm25",\
               "PM10" :"pm10",\
               "SO2"  :"so2"}

input_paths = {value: root_dir+value+"/" for value in species_dir}
input_tails = {value: species_str[value]+root_tail for value in species_dir}

#output_paths = input_paths
output_path = root_dir+"new_anthro_emiss_emissions_netcdf/"


header_end = "NODATA"

var_names = ["agric","domcom","energyprod","indcom","indproc","nature",\
             "offshore","othertrans","roadtrans","solvents","total",\
             "totarea","waste"]


### anthro_emiss specific attributes (make setting these more interactive later?)

netcdf_Att = {}
netcdf_Att['units']    = 'tonnes/km^2/year'

netcdf_sect_map = {}

netcdf_Global_Att = {}



##### point source configuration information

## path to excel spreadsheet with point source information
ps_file = root_dir+'download_data/'+'NAEIPointsSources_2016.xlsx'


## sector mapping for point sources
sector_mapping = {}
sector_mapping['Iron & steel industries']                   = 'indcom'
sector_mapping['Waste collection, treatment & disposal']    = 'waste'
sector_mapping['Other industries']                          = 'indproc'
sector_mapping['Paper, printing & publishing industries']   = 'indproc'
sector_mapping['Other mineral industries']                  = 'indproc'
sector_mapping['Vehicles']                                  = 'roadtrans'
sector_mapping['Oil & gas exploration and production']      = 'offshore'
sector_mapping['Non-ferrous metal industries']              = 'indcom'
sector_mapping['Food, drink & tobacco industry']            = 'indproc'
sector_mapping['Textiles, clothing, leather & footwear']    = 'indproc'
sector_mapping['Chemical industry']                         = 'solvents'
sector_mapping['Other fuel production']                     = 'indproc'
sector_mapping['Major power producers']                     = 'energyprod'
sector_mapping['Electrical engineering']                    = 'indproc'
sector_mapping['Lime']                                      = 'indproc'
sector_mapping['Cement']                                    = 'indproc'
sector_mapping['Public administration']                     = 'domcom'
sector_mapping['Processing & distribution of petroleum products'] = 'offshore'
sector_mapping['Processing & distribution of natural gas']  = 'offshore'
sector_mapping['Agriculture, forestry & fishing']           = 'agric'
sector_mapping['Commercial']                                = 'domcom'
sector_mapping['Water & sewerage']                          = 'waste'
sector_mapping['Minor power producers']                     = 'energyprod'
sector_mapping['Mechanical engineering']                    = 'indproc'
sector_mapping['Miscellaneous']                             = 'indcom'
sector_mapping['Construction']                              = 'indproc'

sector_list = ['energyprod','domcom','indcom','indproc','offshore',\
               'solvents','roadtrans','othertrans','waste','agric','nature']

chem_species = ['Ammonia','Carbon monoxide','Oxides of nitrogen','Sulpher dioxide',\
                'Non-methane VOC','PM10','PM2.5']
chem_mapping = {}
chem_mapping['Non-methane VOC'] = 'NMVOC'
chem_mapping['Sulpher dioxide'] = 'SO2'
chem_mapping['Carbon monoxide'] = 'CO'
chem_mapping['Ammonia'] = 'NH3'
chem_mapping['Oxides of nitrogen'] = 'NOx'
chem_mapping['PM10'] = 'PM10'
chem_mapping['PM2.5'] = 'PM25'



##### functions


def write_to_netcdf(emiss_data,species_list,output_path,netcdf_names):
    """
        simple function for writing individual netcdf files
    """

    for species in species_list:
    
        # copy dataset template
        ds_working = emiss_data[species]

        # write to netcdf file
        outfile = output_path+netcdf_names[species]+'_emiss.nc'
        ds_working.to_netcdf(path=outfile,mode='w')




###### processing code

if __name__=='__main__':

    ## load the emissions data - returning a dictionary containing all datasets
    print('==== Initiating Loading Emission Data ====')
    emiss_data = can.load_emission_data(species_dir,species_str,input_paths,var_names,input_tails,\
                                                        header_end,netcdf_Att,netcdf_sect_map,netcdf_Global_Att,\
                                                        netcdf_names)


    ## point source processing
    print('==== Initiating Point Source Processing ====')
    emiss_data = psa.point_source_processing(emiss_data,ps_file,chem_species,sector_mapping,\
                                                                chem_mapping,species_dir,sector_list)



    ## write emissions data to file
    print('==== Writing emissions to netcdf files ====')
    write_to_netcdf(emiss_data,species_dir,output_path,netcdf_names)




















