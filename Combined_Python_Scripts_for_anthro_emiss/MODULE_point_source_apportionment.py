# -*- coding: utf-8 -*-


#from pathlib import Path
#import os
import numpy as np
import xarray as xr
import pandas as pd



#%% setting the sector mapping for the emissions



def find_nearest(array, value):
    """
    https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx][0]  # we assume that only 1 value is found, this might only be true for 1D arrays?




def extract_emission_grid_info(emiss_data):
    """
    extract the arrays of northing and easting data, as well as a
    template array to use for create chemical temporary data arrays
    """
    
    # just get the information from the first dataset
    for key in emiss_data:
        northing_centre = emiss_data[key].northings.copy(deep=True).to_dataframe()
        easting_centre  = emiss_data[key].eastings.copy(deep=True).to_dataframe()
        break

    easting_grid, northing_grid = np.meshgrid(easting_centre,northing_centre)

    template_grid  = pd.DataFrame(easting_grid)

    return(template_grid,easting_centre,northing_centre)



def create_temporary_chemical_dict(grid_template,chem_list,sector_list):
    """
       creating the temporary chemical storage arrays
    """
    
    
    template_array = grid_template * 0.0

    chem_data = {}
    for chem in chem_list:
        chem_data[chem] = {}
        for sector in sector_list:
            chem_data[chem][sector] = template_array.copy(deep=True)

    return(chem_data)



def load_point_data(ps_file,northing_centre,easting_centre,chem_species):
    """
       loading point data from excel, 
       and add information about the grid position, 
       and select only the chemical data we want
    """

    # open excel spreadsheet, and extract the data sheet
    xlspread = pd.ExcelFile(ps_file)
    point_src = xlspread.parse('Data')

    # get rid of unneeded fields
    psj = point_src.drop(['Year','Operator','Region','Unit','Site','Data Type','PollutantD','SectorID'],axis=1)

    # assuming that any data with North/East grid points = 0 are erroneous, so get rid of them
    #psj_dropped_bad_location = psj[psj.Northing == 0.0]
    psj = psj[psj.Northing != 0.0]

    # scan for the the nearest geographic grid centre for each data point
    psj['Northing_Grid'] = psj['Northing'].apply(lambda x: find_nearest(northing_centre, x))
    psj['Easting_Grid']  = psj['Easting'].apply(lambda x: find_nearest(easting_centre, x))

    # setting indexes, and sorting data
    psj_ne = psj.set_index(['Northing_Grid','Easting_Grid','Sector','PlantID','Pollutant']).sort_index()
    # select only the chemical data that we want (from the list given above)
    psj_ne_chem = psj_ne.loc(axis=0)[:,:,:,:,chem_species]

    return(psj_ne_chem)


def gridding_point_source_data(chem_data,psj_ne_chem,northing_centre,easting_centre,sector_mapping,chem_mapping):
    """
        gridding the point source data onto temporary data grids
    """

    for east_index in psj_ne_chem.index.get_level_values(1).unique():
        print(east_index)
        east_pos = easting_centre[easting_centre.eastings==east_index].index.values.astype(int)[0]
    
        psj_temp_east = psj_ne_chem.loc(axis=0)[:,east_index,:,:,:]
    
        for north_index in psj_temp_east.index.get_level_values(0).unique():
            print('\t\t',north_index)
            north_pos = northing_centre[northing_centre.northings==north_index].index.values.astype(int)[0]
        
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

    return(chem_data)



def copy_data_to_arrays(nc_head,emiss_data,sector_list,chem_data):
    """
    Copying emissions data from temporary data grids into xarray data arrays
    
    During this process we will recalculate the total area sources, as well as the total sources,
    in order that these can be compared to the values in the original files to make sure
    the above processes have not added / removed anything that is should not have.
    """             
    for chem in nc_head:
        print("processing ",chem)
        total   = emiss_data[chem]['total']
        totarea = emiss_data[chem]['totarea']
        total[:,:,:]   = 0.0
        totarea[:,:,:] = 0.0
        for sector in sector_list:
            tdata = emiss_data[chem][sector]
            totarea += tdata
            tdata[0,:,:] += chem_data[chem][sector]            
            tdata[1,:,:] += chem_data[chem][sector]
            total += tdata
    
        #my_new_file = nc_work + chem + tailstr
        #orig_files[chem].to_netcdf(my_new_file,format='NETCDF3_CLASSIC')

    return(emiss_data)


def point_source_processing(emiss_data,ps_file,chem_species,sector_mapping,chem_mapping,chem_list,sector_list):
    """
    Function organising the loading, reading, and mapping of point source data.
    """
    
    (template_grid,easting_centre,northing_centre) = extract_emission_grid_info(emiss_data)
    
    chem_data = create_temporary_chemical_dict(template_grid,chem_list,sector_list)
    
    point_data = load_point_data(ps_file,northing_centre,easting_centre,chem_species)
    
    chem_data = gridding_point_source_data(chem_data,point_data,northing_centre,easting_centre,\
                                                                     sector_mapping,chem_mapping)
    
    emiss_data = copy_data_to_arrays(chem_list,emiss_data,sector_list,chem_data)
    
    
    return(emiss_data)



