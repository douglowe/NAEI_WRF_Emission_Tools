# -*- coding: utf-8 -*-

### CSF3 setup instructions
# module load apps/anaconda3/5.2.0/bin
# pip install --user bng_latlon
#
 


from pathlib import Path
import os
#import shutil
#import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd
#from bng_to_latlon import OSGB36toWGS84


op_root = '/mnt/iusers01/support/mbessdl2/Projects/Emissions_Processing/Emission_Datasets/'
#op_root = '/Users/mbessdl2/work/manchester/EMEP/example_emissions_processing/Emission_Datasets/'

## path to excel spreadsheet with point source information
ps_file = op_root+'NAEI_2016/download_data/'+'NAEIPointsSources_2016.xlsx'


## path to netcdf files containing area emission information
nc_origin = op_root+'NAEI_2016/emissions_netcdf/original_emep/'
nc_work = op_root+'NAEI_2016/emissions_netcdf/'

nc_head = ['ch4', 'voc', 'so2', 'co', 'nh3', 'nox', 'pmco', 'pm25']
tailstr = '_emiss.nc'


#nc_files = [x + tailstr for x in nc_head]
#nc_fpath = [nc_work + x for x in nc_files]
#print(nc_files)

#%% setting the sector mapping for the emissions

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
chem_list = nc_head
chem_mapping = {}
chem_mapping['Non-methane VOC'] = 'voc'
chem_mapping['Sulpher dioxide'] = 'so2'
chem_mapping['Carbon monoxide'] = 'co'
chem_mapping['Ammonia'] = 'nh3'
chem_mapping['Oxides of nitrogen'] = 'nox'
chem_mapping['PM10'] = 'pmco'
chem_mapping['PM2.5'] = 'pm25'


#%% housekeeping, backing up existing output files so we don't accidently delete these 
for my_head in nc_head:
   
    my_file = my_head + tailstr
    
    my_new_file = nc_work + my_file
	# check that we don't have existing copies of the netcdf files, if there are then back them up
    if Path(my_new_file).is_file():
        print('backing up file: '+my_file)
        os.rename(my_new_file,my_new_file+'.bckup')
	
#    my_old_file = nc_origin + my_file
#    # copy the original netcdf files to new location (to avoid any double counting!!!)
#    if Path(my_old_file).is_file():
#        print('copying file: '+my_file)
#        shutil.copy2(my_old_file,my_new_file)
#    else:
#        print('file does not exist: '+my_file)

#%% loading the original netcdf datafiles (we will *NOT* be writing to these)

orig_files = {}

for nhead in nc_head:
    npath = nc_origin + nhead + tailstr
    orig_files[nhead] = xr.open_mfdataset(npath)
    
    


#%% functions for use below

# https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx][0]  # we assume that only 1 value is found, this might only be true for 1D arrays?



#%% set northing, easting mapping

northing_centre = orig_files['ch4'].north.copy(deep=True).to_dataframe()
easting_centre  = orig_files['ch4'].east.copy(deep=True).to_dataframe()

#easting_edges = easting_centre - 500.0
#temp_array = pd.DataFrame([easting_edges.east.values[-1]+1000.0],columns=['east'],index=[easting_edges.shape[0]])
#easting_edges = easting_edges.append(temp_array)

#northing_edges = northing_centre - 500.0
#temp_array = pd.DataFrame([northing_edges.north.values[-1]+1000.0],columns=['north'],index=[northing_edges.shape[0]])
#northing_edges = northing_edges = northing_edges.append(temp_array)

easting_grid, northing_grid = np.meshgrid(easting_centre,northing_centre)

easting_grid  = pd.DataFrame(easting_grid)
northing_grid = pd.DataFrame(northing_grid)


#%% creating temporary chemical storage arrays

template_array = easting_grid * 0.0

chem_data = {}
for chem in chem_list:
    chem_data[chem] = {}
    for sector in sector_list:
        chem_data[chem][sector] = template_array.copy(deep=True)


#%% open excel spreadsheet, and extract the data sheet
xlspread = pd.ExcelFile(ps_file)
point_src = xlspread.parse('Data')


#%% sorting, and filtering point source data

# get rid of unneeded fields
psj = point_src.drop(['Year','Operator','Region','Unit','Site','Data Type','PollutantD','SectorID'],axis=1)

# assuming that any data with North/East grid points = 0 are erroneous, so get rid of them
psj_dropped_bad_location = psj[psj.Northing == 0.0]
psj = psj[psj.Northing != 0.0]

# scan for the the nearest geographic grid centre for each data point
psj['Northing_Grid'] = psj['Northing'].apply(lambda x: find_nearest(northing_centre, x))
psj['Easting_Grid']  = psj['Easting'].apply(lambda x: find_nearest(easting_centre, x))

# setting indexes, and sorting data
psj_ne = psj.set_index(['Northing_Grid','Easting_Grid','Sector','PlantID','Pollutant']).sort_index()
# select only the chemical data that we want (from the list given above)
psj_ne_chem = psj_ne.loc(axis=0)[:,:,:,:,chem_species]


#%% gridding the point source data onto temporary data grids

for east_index in psj_ne_chem.index.get_level_values(1).unique():
    print(east_index)
    east_pos = easting_centre[easting_centre.east==east_index].index.values.astype(int)[0]
    
    psj_temp_east = psj_ne_chem.loc(axis=0)[:,east_index,:,:,:]
    
    for north_index in psj_temp_east.index.get_level_values(0).unique():
        print('\t\t',north_index)
        north_pos = northing_centre[northing_centre.north==north_index].index.values.astype(int)[0]
        
        psj_temp_ne = psj_temp_east.loc(axis=0)[north_index,:,:,:,:]
        #print(psj_ne_chem.loc(axis=0)[north_index,east_index,:])
        
        for sect_index in psj_temp_ne.index.get_level_values(2).unique():
            if(sect_index in sector_mapping):
                psj_temp_sec = psj_temp_ne.loc(axis=0)[:,:,sect_index,:,:]
                
                sector = sector_mapping[sect_index]
        
                for site_id in psj_temp_sec.index.get_level_values(3).unique():
                    psj_temp_site = psj_temp_sec.loc(axis=0)[:,:,:,site_id,:]
                    
                    for poll_index in psj_temp_site.index.get_level_values(4).unique():
                        if(poll_index in chem_mapping):
                            
                            chem = chem_mapping[poll_index] 
                            
                            chem_data[chem][sector].values[north_pos,east_pos] += \
                                        psj_temp_site.loc(axis=0)[:,:,:,:,poll_index].Emission.values[0]
                    
                        else:
                            print('skipping chemical species: ', poll_index)

            else:
                print('skipping sector: ',sect_index)


#%% copying emissions data from temporary data grids into xarray data arrays, and write to new files
#
#  During this process we will recalculate the total area sources, as well as the total sources,
#  in order that these can be compared to the values in the original files to make sure
#  the above processes have not added / removed anything that is should not have.
                
for chem in nc_head:
    print("processing ",chem)
    total   = orig_files[chem]['total'].load()
    totarea = orig_files[chem]['totarea'].load()
    total[:,:,:]   = 0.0
    totarea[:,:,:] = 0.0
    for sector in sector_list:
        tdata = orig_files[chem][sector].load()
        totarea += tdata
        tdata[0,:,:] += chem_data[chem][sector]            
        tdata[1,:,:] += chem_data[chem][sector]
        total += tdata
    
    my_new_file = nc_work + chem + tailstr
    orig_files[chem].to_netcdf(my_new_file,format='NETCDF3_CLASSIC')

        

        
        
        
