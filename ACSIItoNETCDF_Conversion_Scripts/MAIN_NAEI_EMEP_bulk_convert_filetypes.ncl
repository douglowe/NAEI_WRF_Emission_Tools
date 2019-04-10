;
;  Program for converting NAEI ascii data to a netcdf format
;
;
load "./acsii_netcdf_convert_routines.ncl"
load "./model_grid_manipulation_routines.ncl"

begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; Configuration Information ;;;;;;;;;;;;;;;;;;;;;;;
root_dir = "../../NAEI_2016/"

species_dir = (/"CH4","CO","HCl","NH3","NMVOC","NOx","PM01","PM1","PM25","PM10","SO2"/)
species_str = (/"ch4","co","hcl","nh3","voc","nox","pm0_1","pm1","pm2_5","pm10","so2"/)
netcdf_names = (/"ch4","co","hcl","nh3","voc","nox","pm0_1","pm1","pm25","pmco","so2"/)  ;; EMEP conforming


input_paths = root_dir+species_dir+"/"
input_tails = species_str+"16.asc"

output_paths = input_paths





header_end = "NODATA"

var_names = (/"agric","domcom","energyprod","indcom",\
			"indproc","nature","offshore","othertrans",\
			"roadtrans","solvents","total","totarea",\
			"waste"/)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; End of Configuration Information ;;;;;;;;;;;;;;;;

print("== Reading headers of input file ==")

;;; load the header information (from the first file, and assume
;;;              all data files for all species have the same headers!) 
headers = read_header(input_paths(0),var_names(0)+input_tails(0),header_end)


;;; create the "lat" and "lon" variables, which will define the data grid
;;;    We create the lat & lon data arrays here, then modify them inside the procedure 
if(isatt(headers,"ncols"))
	ncols = headers@ncols
else if(isatt(headers,"NCOLS"))
	ncols = headers@NCOLS
else
	print("Number of columns is not defined (need 'ncols' or 'NCOLS' header).")
	exit()
end if
end if
if(isatt(headers,"nrows"))
	nrows = headers@nrows
else if(isatt(headers,"NROWS"))
	nrows = headers@NROWS
else
	print("Number of columns is not defined (need 'ncols' or 'NCOLS' header).")
	exit()
end if
end if
lat_data = new((/nrows,ncols/),double,-9999)
lon_data = new((/nrows,ncols/),double,-9999)
north_data = new(nrows,double,-9999)
east_data  = new(ncols,double,-9999)

print("== Creating the Lat/Lon grid from the header information ==")

calculate_lat_lon_grids(lat_data,lon_data,north_data,east_data,headers)


print("== Loop through emission species ==")
do spci = 0, dimsizes(species_str)-1


	print("=== Creating the output file for species '"+species_str(spci)+"' ===")

	;;; create the netcdf file, and dimensions, global attributes, and latitude & longitude data
	OUTFILE = netcdf_file_create(output_paths(spci),netcdf_names(spci)+"_emiss.nc",\
								lat_data,lon_data,north_data,east_data,headers,var_names,"latitude","longitude")


	;;; process each of the data sets we're loading
	do vari = 0, dimsizes(var_names)-1
		print("==== Processing variable '"+var_names(vari)+"' ====")

		;; test to make sure source file exists, if not then ignore it (but tell us)
		test_file = input_paths(spci)+var_names(vari)+input_tails(spci)
		status = systemfunc(" [ -r "+test_file+" ] && echo 'True' || echo 'False' ")
		if(status .eq. "True")

			;;; load the data table
			data_table := read_table(input_paths(spci),var_names(vari)+input_tails(spci),"float",headers)

			;;; flip the data table vertically (to get correct match with latitudes)
			data_table := data_table(::-1,:)

			;;; save data
			save_variable_data(OUTFILE,data_table,var_names(vari),headers)
		else
			print(test_file+" does not exist - set emissions for this sector to zero")
			data_table := new((/nrows,ncols/),float,-9999)
			data_table(:,:) = 0.0
			;;; save data
			save_variable_data(OUTFILE,data_table,var_names(vari),headers)	
		end if 


	end do


end do

end ; end program
