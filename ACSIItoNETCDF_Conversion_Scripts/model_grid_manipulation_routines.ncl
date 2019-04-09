;
;  Functions for manipulating the model grid information. This includes
;  creating the northings & eastings model grid.
;
;


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./BNGtolatlonconv.ncl"





;
; calculate the latitude and longitude grids for northing and easting grids
;   created from the model data information  
;
undef("calculate_lat_lon_grids")
procedure calculate_lat_lon_grids(lat_data:numeric,lon_data:numeric,north_data:numeric,east_data:numeric,head:logical)
local n_data, e_data, eastings, northings, options, lat_lon, ncols, nrows
begin

; pull out the information we need from the headers
if(isatt(head,"ncols"))
	ncols = head@ncols
else if(isatt(head,"NCOLS"))
	ncols = head@NCOLS
else
	print("Number of columns is not defined (need 'ncols' or 'NCOLS' header). Stopping in calculate_lat_lon_grids. ")
	exit()
end if
end if

if(isatt(head,"nrows"))
	nrows = head@nrows
else if(isatt(head,"NROWS"))
	nrows = head@NROWS
else
	print("Number of columns is not defined (need 'ncols' or 'NCOLS' header). Stopping in calculate_lat_lon_grids. ")
	exit()
end if
end if

if(isatt(head,"xllcorner"))
	xllcorner = head@xllcorner
else if(isatt(head,"XLLCORNER"))
	xllcorner = head@XLLCORNER
else
	print("X Lower-Left Corner is not defined (need 'xllcorner' or 'XLLCORNER' header). Stopping in calculate_lat_lon_grids. ")
	exit()
end if
end if

if(isatt(head,"yllcorner"))
	yllcorner = head@yllcorner
else if(isatt(head,"YLLCORNER"))
	yllcorner = head@YLLCORNER
else
	print("Y Lower-Left Corner is not defined (need 'yllcorner' or 'YLLCORNER' header). Stopping in calculate_lat_lon_grids. ")
	exit()
end if
end if

if(isatt(head,"cellsize"))
	cellsize = head@cellsize
else if(isatt(head,"CELLSIZE"))
	cellsize = head@CELLSIZE
else
	print("Grid cell size is not defined (need 'cellsize' or 'CELLSIZE' header). Stopping in calculate_lat_lon_grids. ")
	exit()
end if
end if



; set up the northing & easting data arrays
n_data = new((/nrows,ncols/),float,-9999)
e_data = new((/nrows,ncols/),float,-9999)

; create 1D easting and westing arrays
eastings  = fspan(xllcorner,(xllcorner + (ncols-1)*cellsize),ncols)
northings = fspan(yllcorner,(yllcorner + (nrows-1)*cellsize),nrows)

; switch from lower-left corner to middle of gridcell
eastings  = eastings  + cellsize/2.0
northings = northings + cellsize/2.0

; copy to 1D arrays for later output
east_data = (/eastings/)
north_data = (/northings/)

; copy the 1D arrays into the 2D arrays
do ii = 0,nrows-1
	e_data(ii,:) = (/eastings/)
end do

do ii = 0,ncols-1
	n_data(:,ii) = (/northings/)
end do

; set up the options for our Lat/Lon conversion 
options = True
options@ellipsoid = "airy1830_ellipsoid"
options@grid = "bng_grid"

lat_lon = utm_latlon(todouble(n_data),todouble(e_data),options)

; copy the lat / lon data to the arrays for output
lat_data = (/lat_lon(0,:,:)/)
lon_data = (/lat_lon(1,:,:)/)


end

