; Routines for converting between British National Grid coordinates and latitude and longitude
;
; Calculations are from Annexe C of 'A guide to coordinate systems in Great Britain' by
; Ordnance Survey.
; https://en.wikipedia.org/wiki/Ordnance_Survey_National_Grid
; http://www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf
;
; Note that the datum used is OSGB 1936. A transformation will need to be applied if
; working with another datum.
;

;
;  NCL Note - need NCL version 6.3.0 or later
;        (for get_unique_values function)
;


;
;  required constants
;
PI =  3.141592653589793d0

;
; constants for Airy 1830 ellipsoid
;
airy1830_ellipsoid = True
airy1830_ellipsoid@a = 6377563.396d0  ; semi-major axis a (m)
airy1830_ellipsoid@b = 6356256.909d0  ; semi-minor axis b (m)
 
;
; constants for British National Grid
;
bng_grid = True
bng_grid@N0 = -100000.0d0		; northing coordinates for true origin (m) 
bng_grid@E0 = 400000.0d0		; easting coordinates for true origin (m)
bng_grid@F0 = 0.9996012717d0	; scale factor on central meridian
bng_grid@Lat_Ref = 49.0d0		; Latitude of true origin (degrees)
bng_grid@Lon_Ref = -2.0d0		; Longitude of true origin (degrees)


;
; Function for calculating "M" (equation C3) given set values of 
;    b, F0, ns, phi and phi0
;
undef("M_calculation")
function M_calculation(b:double,F0:double,ns:double,phi:double,phi0:double)
local M1, M2, M3, M4
begin

	M1 = (1.0d0+ns+((5.0d0/4.0d0)*(ns^2.0d0))+((5.0d0/4.0d0)*(ns^3.0d0))) \
						* (phi-phi0)
	M2 = ((3.0d0*ns)+(3.0d0*(ns^2.0d0))+((21.0d0/8.0d0)*(ns^3.0d0))) \
						* (sin(phi-phi0)) \
						* (cos(phi+phi0))
	M3 = (((15.0d0/8.0d0)*(ns^2.0d0))+((15.0d0/8.0d0)*(ns^3.0d0))) \
						* (sin(2.0d0*(phi-phi0))) \
						* (cos(2.0d0*(phi+phi0)))
	M4 = ((35.0d0/24.0d0)*(ns^3.0d0)) \
						* (sin(3.0d0*(phi-phi0))) \
						* (cos(3.0d0*(phi+phi0)))

	M = b * F0 * (M1-M2+M3-M4)
	
	return(M)

end




;
;  Function for calculating phi_prime for unique Northing values
;
;  This returns an array the same size as the Northing array, containing
;     these values. This array will be of the type "double", otherwise
;     too much information will be lost. Likewise, all the calculations
;     are made using double precision.
;
undef("phi_prime_iteration")
function phi_prime_iteration(norths:double,ellipsoid:logical,grid:logical)
local phi0, ns, phi_prime_array, north_unique, \
		phi_prime, cond, M1, M2, M3, M4, M
begin

	; calculate phi0 for grid projection
	phi0 = todouble(grid@Lat_Ref) * (PI/180.0d0)
	; calculate ns for the ellipsoid
	ns = (todouble(ellipsoid@a) - todouble(ellipsoid@b)) / (todouble(ellipsoid@a) + todouble(ellipsoid@b))

	; create phi_prime data array
	phi_prime_array = new(dimsizes(norths),double,-9999.)

	; find the unique northings values
	north_unique = get_unique_values(norths)
	
	; loop through the unique values
	do ii = 0, dimsizes(north_unique)-1
		
		; initialise phi_prime to phi0 (needed for first step)
		phi_prime = phi0
		
		; initialise the exit condition variables
		cond = 1.0d0	; (has to be >0.01e-3 for the first step)
		M = 0.0d0		; (has to be 0.0 for the first step)

		; iterating through the calculations for Equation C7
		do while (cond .ge. 0.001d-3)
			; calculate phi_prime
			phi_prime = ( ( (/north_unique(ii)/) - grid@N0 - M ) / ( todouble(ellipsoid@a) * todouble(grid@F0) ) ) \
																	+ phi_prime

			; calculating Equation C3
			M = M_calculation(todouble(ellipsoid@b),todouble(grid@F0),ns,phi_prime,phi0)
		
			; calculate the exit condition for the next loop
			cond = north_unique(ii) - todouble(grid@N0) - M
			
		end do
				
		; record the phi_prime value for all locations matching the northing value
		phi_prime_array = where(norths.eq.north_unique(ii),phi_prime,phi_prime_array)

	end do
	
	
	return(phi_prime_array)

end









;
;  Main function for organising the calculation of the new Lat / Lon values
;
;  To use this function we need to pass this information:
; 
;  norths --- 1D or 2D array, containing northings data
;  easts  --- 1D or 2D array, containing eastings data
;  options --- logical for selecting the ellipsoid and grid to use
; 
;  ellipsoid --- logical, with attributes for the global ellipsoid:
;           a -- semi-major axis a (m)
;           b -- semi-minor axis b (m)
;  grid --- logical, with attributes for the transverse mercator projection
;           N0 -- northing coordinates for true origin (m) 
;           E0 -- easting coordinates for true origin (m) 
;           F0 -- scale factor on central meridian
;           Lat_Ref -- Latitude of true origin (degrees)
;           Lon_Ref -- Longitude of true origin (degrees)
;
;  This routine returns:
;   lat_lon --- 2D or 3D array, containing the lat and lon data
;               	The first index will be for lat/lon 
;					(lat will be value 0, lon will be value 1) 
;
undef("utm_latlon")
function utm_latlon(norths:double,easts:double,options:logical)
local nw_size, nw_ndim, lat_lon, lat0, lon0, ii, jj, \
		lambda0, phi_prime, e2, v, rho, eta2, \
		VII, VIII, IX, X, XI, XII, XIIA, phi, lambda, \
		grid, ellipsoid
begin

	; set the ellipsoid and grid information
	if(options)
		if(isatt(options,"ellipsoid"))
			if(options@ellipsoid .eq. "airy1830_ellipsoid")
				ellipsoid = airy1830_ellipsoid		
			else
				print("ellipsoid option '"+options@ellipsoid+"' not recognised, using airy1830_ellipsoid instead")			
				ellipsoid = airy1830_ellipsoid		
			end if
		else
			print("ellipsoid option not set, using airy1830_ellipsoid default option")			
			ellipsoid = airy1830_ellipsoid		
		end if

		if(isatt(options,"grid"))
			if(options@grid .eq. "bng_grid")
				grid = bng_grid		
			else
				print("gird option '"+options@grid+"' not recognised, using bng_grid instead")			
				grid = bng_grid
			end if
		else
			print("grid option not set, using bng_grid default option")			
			grid = bng_grid		
		end if

	else ; if not set then use defaults
		ellipsoid = airy1830_ellipsoid
		grid      = bng_grid
	end if

	; get the size of the original norths / easts data arrays (assume that they are the same size!)
	nw_size = dimsizes(norths)
	; get the number of dimensions
	nw_ndim = dimsizes(nw_size)
	
	; create the new lat_lon array, depending on data size
	if(nw_ndim .eq. 1)
		lat_lon = new((/2,nw_size(0)/),double,-9999.)
	else if(nw_ndim .eq. 2)
		lat_lon = new((/2,nw_size(0),nw_size(1)/),double,-9999.)	
	else
		print("function utm_latlon can't cope with input arrays with "\
				+nw_ndim+" dimensions, please use 1 or 2 dimensions instead")
		exit()
	end if
	end if
		
	
	; calculate phi_prime for every unique northing value (Equation C3)
	phi_prime = phi_prime_iteration(norths,ellipsoid,grid)
	
	; calculate e2 for the ellipsoid
	e2 = ((todouble(ellipsoid@a)^2.0d0)-(todouble(ellipsoid@b)^2.0d0)) / (todouble(ellipsoid@a)^2.0d0)
	
	; calculate the values for v, rho and eta2, Equations C2
	v = todouble(ellipsoid@a) * todouble(grid@F0) * \
		( (1.0d0-(e2*(sin(phi_prime)^2.0d0)))^(-0.5d0) )
	rho = todouble(ellipsoid@a) * todouble(grid@F0) * (1.0d0-e2) * \
		( (1.0d0-(e2*(sin(phi_prime)^2.0d0)))^(-1.5d0) )
	eta2 = (v/rho)-1.0d0

	; calculate factors for Equations C8 and C9
	VII  = (tan(phi_prime)) / (2.0d0*rho*v)

	VIII = ( (tan(phi_prime)) / (24.0d0*rho*(v^3.0d0)) ) * \
				( 5.0d0 + (3.0d0*(tan(phi_prime)^2.0d0)) + eta2 - (9.0d0*(tan(phi_prime)^2.0d0)*eta2) )

	IX   = ( (tan(phi_prime)) / (720.0d0*rho*(v^5.0d0)) ) * \
				( 61.0d0 + (90.0d0*(tan(phi_prime)^2.0d0)) + (45.0d0*(tan(phi_prime)^4.0d0)) )

	X    = ( 1.0d0 / (cos(phi_prime)) ) / v

	XI   = ( (1.0d0/(cos(phi_prime))) / (6.0d0*(v^3.0d0)) ) * \
				( (v/rho) + (2.0d0*(tan(phi_prime)^2.0d0)) )

	XII  = ( (1.0d0/(cos(phi_prime))) / (120.0d0*(v^5.0d0)) ) * \
			( 5.0d0 + (28.0d0*(tan(phi_prime)^2.0d0)) + (24.0d0*(tan(phi_prime)^4.0d0)) )

	XIIA = ( (1.0d0/(cos(phi_prime))) / (5040.0d0*(v^7.0d0)) ) * \
			( 61.0d0 + (662.0d0*(tan(phi_prime)^2.0d0)) + \
					(1320.0d0*(tan(phi_prime)^4.0d0)) + \
					(720.0d0*(tan(phi_prime)^6.0d0)) )

	; calculate lambda0 for grid projection
	lambda0 = todouble(grid@Lon_Ref) * (PI/180.0d0)

	
	; calculate phi and lambda (Equations C8 and C9)
	phi = phi_prime - (VII*((easts-todouble(grid@E0))^2.0d0)) + \
			(VIII*((easts-todouble(grid@E0))^4.0d0)) - (IX*((easts-todouble(grid@E0))^6.0d0))
	
	lambda = lambda0 + (X*(easts-todouble(grid@E0))) - (XI*((easts-todouble(grid@E0))^3.0d0)) + \
						(XII*((easts-todouble(grid@E0))^5.0d0)) - (XIIA*((easts-todouble(grid@E0))^7.0d0))

	
	; create temporary variables for lat & lon data
	lat0 = new(1,double,-9999.)
	lon0 = new(1,double,-9999.)
	
	; loop through the data arrays, and call the lat lon calculation procedure
	if(nw_ndim .eq. 1)
		; latitude calculation
		lat_lon(0,:) = phi * (180.0d0/PI)
		; longitude calculation
		lat_lon(1,:) = lambda * (180.0d0/PI)
	else if(nw_ndim .eq. 2)
		; latitude calculation
		lat_lon(0,:,:) = phi * (180.0d0/PI)
		; longitude calculation
		lat_lon(1,:,:) = lambda * (180.0d0/PI)
	end if
	end if

	
	return(lat_lon)

end






;
;  Main function for organising the calculation of the new easting and northing values
;
;  To use this function we need to pass this information:
; 
;  lats   --- 1D or 2D array, containing latitude data
;  longs  --- 1D or 2D array, containing longitude data
;  options --- logical for selecting the ellipsoid and grid to use
; 
;  ellipsoid --- logical, with attributes for the global ellipsoid:
;           a -- semi-major axis a (m)
;           b -- semi-minor axis b (m)
;  grid --- logical, with attributes for the transverse mercator projection
;           N0 -- northing coordinates for true origin (m) 
;           E0 -- easting coordinates for true origin (m) 
;           F0 -- scale factor on central meridian
;           Lat_Ref -- Latitude of true origin (degrees)
;           Lon_Ref -- Longitude of true origin (degrees)
;
;  This routine returns:
;   north_east --- 2D or 3D array, containing the northing and easting data
;               	The first index will be for northing/easting 
;					(northing will be value 0, easting will be value 1) 
;
undef("utm_northeast")
function utm_northeast(lats:numeric,lons:numeric,options:logical)
local nw_size, nw_ndim, lat_lon, lat0, lon0, ii, jj, \
		lambda0, phi_prime, e2, v, rho, eta2, \
		sec_I, sec_II, sec_III, sec_IIIa, sec_IV, sec_V, sec_VI, \
		phi, lambda, grid, ellipsoid
begin

	; set the ellipsoid and grid information
	if(options)
		if(isatt(options,"ellipsoid"))
			if(options@ellipsoid .eq. "airy1830_ellipsoid")
				ellipsoid = airy1830_ellipsoid		
			else
				print("ellipsoid option '"+options@ellipsoid+"' not recognised, using airy1830_ellipsoid instead")			
				ellipsoid = airy1830_ellipsoid		
			end if
		else
			print("ellipsoid option not set, using airy1830_ellipsoid default option")			
			ellipsoid = airy1830_ellipsoid		
		end if

		if(isatt(options,"grid"))
			if(options@grid .eq. "bng_grid")
				grid = bng_grid		
			else
				print("gird option '"+options@grid+"' not recognised, using bng_grid instead")			
				grid = bng_grid
			end if
		else
			print("grid option not set, using bng_grid default option")			
			grid = bng_grid		
		end if

	else ; if not set then use defaults
		ellipsoid = airy1830_ellipsoid
		grid      = bng_grid
	end if

	; get the size of the original norths / easts data arrays (assume that they are the same size!)
	nw_size = dimsizes(lats)
	; get the number of dimensions
	nw_ndim = dimsizes(nw_size)
	
	; create the new lat_lon array, depending on data size
	if(nw_ndim .eq. 1)
		north_east = new((/2,nw_size(0)/),double,-9999.)
	else if(nw_ndim .eq. 2)
		north_east = new((/2,nw_size(0),nw_size(1)/),double,-9999.)	
	else
		print("function utm_northeast can't cope with input arrays with "\
				+nw_ndim+" dimensions, please use 1 or 2 dimensions instead")
		exit()
	end if
	end if
		
	
	
	; calculate lambda0 for grid projection
	lambda0 = todouble(grid@Lon_Ref) * (PI/180.0d0)
	; calculate lambda for all longitudes
	lambda = lons * (PI/180.0d0)

	; calculate phi0 for grid projection
	phi0 = todouble(grid@Lat_Ref) * (PI/180.0d0)
	; calculate phi for all latitudes
	phi = lats * (PI/180.0d0)

	; calculate variables in equations C1
	ns = ( todouble(ellipsoid@a) - todouble(ellipsoid@b) ) / ( todouble(ellipsoid@a) + todouble(ellipsoid@b) )
	
	
	; calculate M factor for every latitude value
	M_fac = M_calculation(todouble(ellipsoid@b),todouble(grid@F0),ns,phi,phi0)
	
	
	; calculate e2 for the ellipsoid
	e2 = ((todouble(ellipsoid@a)^2.0d0)-(todouble(ellipsoid@b)^2.0d0)) / (todouble(ellipsoid@a)^2.0d0)
	
	; calculate the values for v, rho and eta2, Equations C2
	v = todouble(ellipsoid@a) * todouble(grid@F0) * \
		( (1.0d0-(e2*(sin(phi)^2.0d0)))^(-0.5d0) )
	rho = todouble(ellipsoid@a) * todouble(grid@F0) * (1.0d0-e2) * \
		( (1.0d0-(e2*(sin(phi)^2.0d0)))^(-1.5d0) )
	eta2 = (v/rho)-1.0d0

	; calculate factors for Equations C4
	sec_I    = M_fac + todouble(grid@N0)
	
	sec_II   = (v/2.0d0) * sin(phi) * cos(phi)

	sec_III  = (v/24.0d0) * sin(phi) * (cos(phi)^3.0d0) * \
				( 5.0d0 - tan(phi)^2.0d0 + 9.0d0*eta2 )

	sec_IIIa = (v/720.0d0) * sin(phi) * (cos(phi)^5.0d0) * \
				( 61.0d0 - 58.0d0 * tan(phi)^2.0d0 + tan(phi)^4.0d0 )

	sec_IV   = v * cos(phi)
	
	sec_V    = (v/6.0d0) * (cos(phi)^3.0d0) * \
				( v/rho - (tan(phi)^2.0d0) )

	sec_VI   = (v/120.0d0) * (cos(phi)^5.0d0) * \
				( 5.0d0 - 18.0d0 * tan(phi)^2.0d0 + tan(phi)^4.0d0 + \
						14.0d0*eta2 - 58.0d0*(tan(phi)^2.0d0)*eta2 )


	; calculate Northings and Eastings (Equations C4 and C5)
	NN = sec_I + \
			sec_II*(lambda-lambda0)^2.0d0 + \
			sec_III*(lambda-lambda0)^4.0d0 + \
			sec_IIIa*(lambda-lambda0)^6.0d0
	
	EE = grid@E0 + \
			sec_IV*(lambda-lambda0) + \
			sec_V*(lambda-lambda0)^3.0d0 + \
			sec_VI*(lambda-lambda0)^5.0d0
	
	
	; loop through the data arrays, and call the lat lon calculation procedure
	if(nw_ndim .eq. 1)
		; latitude calculation
		north_east(0,:) = NN
		; longitude calculation
		north_east(1,:) = EE
	else if(nw_ndim .eq. 2)
		; latitude calculation
		north_east(0,:,:) = NN
		; longitude calculation
		north_east(1,:,:) = EE
	end if
	end if

	
	return(north_east)

end






