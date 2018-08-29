;
;  Program for converting NAEI ascii data to a netcdf format
;
;
load "./acsii_netcdf_convert_routines.ncl"
load "./model_grid_manipulation_routines.ncl"

begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; Configuration Information ;;;;;;;;;;;;;;;;;;;;;;;
input_path = "../nmvoc_data/"
input_tail = "voc14.asc"

output_path = "./"
netcdf_name = "nmvoc"


header_end = "NODATA"

var_names = (/"agric","domcom","energyprod","indcom",\
			"indproc","nature","offshore","othertrans",\
			"roadtrans","solvents","total","totarea",\
			"waste"/)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; End of Configuration Information ;;;;;;;;;;;;;;;;

print("== Reading headers of input file ==")

;;; load the header information (from the first file, and assume
;;;              all data files have the same headers! 
headers = read_header(input_path,var_names(0)+input_tail,header_end)


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

print("== Creating the Lat/Lon grid from the header information ==")

calculate_lat_lon_grids(lat_data,lon_data,headers)


print("== Creating the output file ==")

;;; create the netcdf file, and dimensions, global attributes, and latitude & longitude data
OUTFILE = netcdf_file_create(output_path,netcdf_name+"_emiss.nc",lat_data,lon_data,headers,var_names)


;;; process each of the data sets we're loading
do vari = 0, dimsizes(var_names)-1
	print("=== Processing variable '"+var_names(vari)+"' ===")

	;;; load the data table
	data_table = read_table(input_path,var_names(vari)+input_tail,"float",headers)

	;;; flip the data table vertically (to get correct match with latitudes)
	data_table = data_table(::-1,:)

	;;; create variables in netcdf file, and save data
	save_variable_data(OUTFILE,data_table,var_names(vari),headers)


end do


end ; end program