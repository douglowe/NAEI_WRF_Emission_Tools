load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;
;  Function for extracting the header data from an ASCII file, and
;     also to count how deep the header is.
;
;  This function requires a string value to be passed, which defines
;     the string starts the final line of the header.
;
undef("read_header")
function read_header(filepath:string,filename:string,headend:string)
local head_strings, head_len, headers, ii
begin
	; extract the header strings
	;head_strings = readAsciiHead("../nmvoc_data/indcomvoc14.asc","NODATA")
	head_strings = readAsciiHead(filepath+filename,headend)
	head_len = dimsizes(head_strings)


	headers = True
	headers@header_length = head_len

	; loop through header, splitting by white space, and saving the 
	;   header information as attributes.
	; Note: we assume that each header line has two values, and that the 2nd is numeric
	;
	do ii = 0,head_len-1
		strs = str_split(head_strings(ii)," ")
		headers@$strs(0)$ = stringtointeger(strs(1))
	end do

	return(headers)

end

;
;  Function for extracting table of data from ASCII file, and then add the 
;     header information to the attributes, and set the missing value number
;
undef("read_table")
function read_table(filepath:string,filename:string,datatype:string,head:logical)
local data_table, ncols, head_length
begin
	; check that we are defining the number of columns to be read
	if(isatt(head,"ncols"))
		ncols = head@ncols
	else if(isatt(head,"NCOLS"))
		ncols = head@NCOLS
	else
		print("Number of columns to be read is not defined (need 'ncols' or 'NCOLS' header), halting program in 'read_table'.")
		exit()
	end if
	end if
	
	; look for header length (or assume this)
	if(isatt(head,"header_length"))
		head_length = head@header_length
	else
		print("Header length not defined (need 'header_length' attribute), halting program in 'read_table'.")
		exit()
	end if
	
	; read data table
	;data_table = readAsciiTable("../nmvoc_data/indcomvoc14.asc",headers@ncols,"float",headers@header_length)
	data_table = readAsciiTable(filepath+filename,ncols,datatype,head_length)

	; check if we have defined a missing data value, if so then replace it
	if(isatt(head,"NODATA_value"))
		data_table@_FillValue = head@NODATA_value
	else
		print("No data value not defined (need 'NODATA_value' header), using default instead")
	end if

	; replace missing values with zeros
	data_table = where(.not.ismissing(data_table),data_table,0.0)

	return(data_table)

end


undef("netcdf_file_create")
function netcdf_file_create(filepath:string,filename:string,lat_data:numeric,lon_data:numeric,\
							north_data:numeric,east_data:numeric,\
							head:logical,var_names:string,lat_name:string, lon_name:string)
local OUTFILE, dimNames, dimSizes, dimUnlim, globalatt, Time_Var, Data_Var, vari, \
		ncols, nrows
begin

	; check if the file exists already, if so then back the old file up
	if(fileexists(filepath+filename))
		system("mv "+filepath+filename+" "+filepath+filename+".backup")
	end if

	; create output file
	OUTFILE = addfile(filepath+filename,"c")
	
	; pull out column and row information
	if(isatt(head,"ncols"))
		ncols = head@ncols
	else if(isatt(head,"NCOLS"))
		ncols = head@NCOLS
	else
		print("Number of columns is not defined (need 'ncols' or 'NCOLS' header), halting program in 'netcdf_file_create'.")
		exit()
	end if
	end if

	if(isatt(head,"nrows"))
		nrows = head@nrows
	else if(isatt(head,"NROWS"))
		nrows = head@NROWS
	else
		print("Number of rows is not defined (need 'nrows' or 'NROWS' header), halting program in 'netcdf_file_create'.")
		exit()
	end if
	end if

	
	;;;; define the dimensions
	;; time -> use 2 (so anthro_emiss will interpolate okay)
	;; lat  -> use number of rows
	;; lon  -> use number of columns
	dimNames = (/"time","lat","lon"/)
	dimSizes = (/2,nrows,ncols/)
	dimUnlim = (/False,False,False/) 
	
	;;;; define the global attributes
	globalAtt = True
	globalAtt@info = "Test file"
	
	;;;; define the time variable
	Date_Var = new(2,integer,-9999)
	Date_Var!0 = "time"
	Date_Var(0) = 19000101
	Date_Var(1) = 21000101
	Date_Var@units = "YYYYMMDD"
	Date_Var@long_name = "Date"

	;;;; name the dimensions for lat and lon data
	lat_data!0 = "lat"
	lon_data!0 = "lat"
	lat_data!1 = "lon"
	lon_data!1 = "lon"
	lat_data@units = "degrees_north"
	lon_data@units = "degrees_east"
	lat_data@long_name = "Latitude"
	lon_data@long_name = "Longitude"
	;;;; name the dimensions for northing and easting data
	north_data!0 = "lat"
	east_data!0  = "lon"
	north_data@units = "metres north from origin"
	east_data@units  = "metres east from origin"
	north_data@long_name = "Northings"
	east_data@long_name  = "Eastings"
	
	;;;; create a dummy data variable
	Data_Var = new((/2,nrows,ncols/),float,-9999)
	Data_Var!0 = "time"
	Data_Var!1 = "lat"
	Data_Var!2 = "lon"
	Data_Var@units = "tonnes/km^2/year"
	
	;;;; define the file settings
	setfileoption(OUTFILE,"DefineMode",True)
	
		; set file attributes
		filedimdef(OUTFILE,dimNames,dimSizes,dimUnlim)
		fileattdef(OUTFILE,globalAtt)
	
		; create the lat and lon data arrays
		filevardef(OUTFILE, lat_name, typeof(lat_data), getvardims(lat_data))                          
		filevardef(OUTFILE, lon_name, typeof(lon_data), getvardims(lon_data))                          
		filevarattdef(OUTFILE, lat_name, lat_data)
		filevarattdef(OUTFILE, lon_name, lon_data)

		; create the northing and easting data arrays
		filevardef(OUTFILE, "north", typeof(north_data), getvardims(north_data))
		filevardef(OUTFILE, "east", typeof(east_data), getvardims(east_data))
		filevarattdef(OUTFILE, "north", north_data)
		filevarattdef(OUTFILE, "east", east_data)
	
		; create the time data array
		filevardef(OUTFILE, "date", typeof(Date_Var), getvardims(Date_Var))
		filevarattdef(OUTFILE, "date", Date_Var)		
		
		; create the data arrays
		do vari = 0, dimsizes(var_names)-1
			filevardef(OUTFILE, var_names(vari), typeof(Data_Var), getvardims(Data_Var))                          
			filevarattdef(OUTFILE, var_names(vari), Data_Var)
		end do
		
	setfileoption(OUTFILE,"DefineMode",False)

	;;;; copy time and position data over to the file
	OUTFILE->$lat_name$ = (/lat_data/)
	OUTFILE->$lon_name$ = (/lon_data/)
	OUTFILE->date = (/Date_Var/)
	OUTFILE->north = (/north_data/)
	OUTFILE->east  = (/east_data/)


	; return output file handle
	return(OUTFILE)


end




undef("save_variable_data")
procedure save_variable_data(OUTFILE:file,data_table:numeric,var_name:string,head:logical)
begin

	;;; copy data to the new variable
	OUTFILE->$var_name$(0,:,:) = (/data_table/)
	OUTFILE->$var_name$(1,:,:) = (/data_table/)

end



