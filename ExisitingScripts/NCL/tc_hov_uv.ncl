; Script to produce Hovmoller plots of three-dimensional windspeed

; Run using

; ncl dat=\"02T12\" opt=\"x11\" ens0=\"em11\" dist=3.0 ts=42 tf=86 mlev0=13 rmw=1 calc=1
; tclr=2 rclr=2 wclr=1 lay=1 mlev1=24 ar=1.5 nr=31 ar0=1.5 cn0=\"geo_sm\" tc_hov_uv.ncl

; 'dat'    = initialisation time string ("02T12", "03T00", ...)
; 'ens0'   = ensemble simulation (em00-em11)
; 'dist'   = size of box following storm (degrees)
; 'opt'    = output file format ("pdf" or "x11")
; 'mlev'   = model level for plotting PV and vertical velocity (see L106-108)
; 'rmw'	   = read in RMW data on model (1) or pressure (0) levels
; 'calc'   = centre calculated on each level (0), calculated offline (1) or using 980 m vort (2)
; 'tclr'   = colour scheme for plotting tangential wind (0-2)
; 'rclr'   = colour scheme for plotting radial wind (0-2)
; 'wclr'   = colour for plotting vertical velocity (0-2)
; 'lay'    = calculate layer average (1) or single level (0)
; 'mlev2'  = upper model level for calculating layer average (only valid when lay=1)
; 'ar'	   = distance in degrees of outer radius (1.0, 1.5, etc; default is 1.5)
; 'nr'	   = number of radial circles between r = 0 and r = 'ar' (21, 41, etc; default is 31)
; 'ar0'	   = outer radius (degrees) for final plot [must be <= 'ar']
; 'cn0'    = method of storm tracking ["slp", "geo", "geo_sm", "vort"] - default is "geo_sm"

; 'mlev'  = 13 and 'mlev2' = 24 give an approximate 1-3 km layer average

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'

dat         = "02T12"
opt         = "pdf"
ens0        = "em11"
dist        = 3.0
ts          = 42
tf          = 86
mlev0       = 13
rmw         = 1
calc        = 1
tclr        = 2
rclr        = 2
wclr        = 1
lay         = 1
mlev1       = 24
ar          = 1.5
nr          = 31
ar0         = 1.5
cn0         = "geo_sm"

; Load main NCL libraries ($NCARG_ROOT = /apps/applications/ncl/6.5.0/1/default)
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/nfs/a37/earshar/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

; Also load additional user-written functions
load "/nfs/a37/earshar/ncl_func/st_rm.ncl"
load "/nfs/a37/earshar/ncl_func/nc_times.ncl"
load "/nfs/a37/earshar/ncl_func/nc_levs.ncl"

begin

;=================
; Read in files
;=================

; Find 'pb' stream data to read in and plot
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_c  = numINFO_c[0]            ; Number of times in 'pc' stream
  time_c      = numINFO_c[1]            ; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream file

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, etc)
; 'ens0'    = ensemble simulation (em00-em11)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 0)
  u_cyc     = storm_rel[0]             ; Zonal wind
  v_cyc     = storm_rel[1]             ; Meridional wind
  vel_cyc   = storm_rel[2]             ; Vector wind
  numTIMES  = storm_rel[3]             ; Number of times in file
  lat_arr   = storm_rel[4]             ; Latitude array subset (following storm)
  lon_arr   = storm_rel[5]             ; Longitude array subset (following storm)
  centre    = storm_rel[6]             ; Storm track (position) array

; Create arrays for data and output information
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;=======================================
; Choose which storm track to read in
;=======================================

  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "vort") then
   cn = 1
  elseif (cn0 .eq. "geo") then
   cn = 2
  elseif (cn0 .eq. "geo_sm") then
   cn = 3
  end if

;==========================================================================
; Create date/time string arrays for all times in file (use 'pc' stream)
;==========================================================================

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date = cd_calendar(time_c, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59)
  do ct = 0, numINPUT_c-1
   if (minute(ct).gt.30) then
     hour(ct) = hour(ct)+1
   end if
  end do
  delete(ct)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_c,string)
  time_str  = new(numINPUT_c,string)
  time_arr  = new(numINPUT_c,string)
  title_arr = new(numINPUT_c,string)

  date_str = sprinti("%0.2iUTC ", hour) + \
             sprinti("%0.2i ", day) \
           + month_abbr(month)
  time_str = sprinti("%0.2i UTC ", hour) + \
             sprinti("%0.2i ", day) \
           + month_abbr(month)
  time_arr = sprinti("%0.2i", day) + \
             month_abbr(month) + \
             "_" + sprinti("%0.2iZ", hour)

;========================================
; Define cylindrical coordinate arrays
;========================================

; ar = distance (degrees) of outer radial circle (1.0)
; nr = total number of radial circles (21)

  thetas     = new(73,float)
  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define arrays for later use in tangential velocity cross-section plots
  vtan_azi   = new((/dimsizes(radii)/),float) ; Tangential wind
  vrad_azi   = new((/dimsizes(radii)/),float) ; Radial wind
  spd_azi    = new((/dimsizes(radii)/),float) ; Vector velocity
  pvrt_azi   = new((/dimsizes(radii)/),float) ; Potential vorticity
  vvel_azi   = new((/dimsizes(radii)/),float) ; Vertical velocity

  nTIMES     = (tf - ts) + 1		      ; Number of times in user-defined subset

; Define 2-dimensional array to hold all values from 'vtan_azi' at ALL times
  vtan_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vtan_plot0!0 = "rad"
  vtan_plot0!1 = "time"

  vrad_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vrad_plot0!0 = "rad"
  vrad_plot0!1 = "time"

  pvrt_plot0   = new((/dimsizes(radii),numTIMES/),float)
  pvrt_plot0!0 = "rad"
  pvrt_plot0!1 = "time"

  vvel_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vvel_plot0!0 = "rad"
  vvel_plot0!1 = "time"

; Add metadata for azimuthally averaged variables
  vtan_azi!0           = "rad"
  vtan_azi&rad         = radii
  vtan_azi@description = "Azimuthally averaged tangential velocity"
  vtan_azi@units       = "m s~S~-1~N~"

  vrad_azi!0           = "rad"
  vrad_azi&rad         = radii
  vrad_azi@description = "Azimuthally averaged radial velocity"
  vrad_azi@units       = "m s~S~-1~N~"

  spd_azi!0            = "rad"
  spd_azi&rad          = radii
  spd_azi@description  = "Azimuthally averaged vector velocity"
  spd_azi@units        = "m s~S~-1~N~"

  pvrt_azi!0           = "rad"
  pvrt_azi&rad         = radii
  pvrt_azi@description = "Azimuthally averaged PV"
  pvrt_azi@units       = "10~S~-6~N~K kg~S~-1~N~ m~S~2~N~ s~S~-1~N~"

  vvel_azi!0            = "rad"
  vvel_azi&rad          = radii
  vvel_azi@description  = "Azimuthally averaged vertical velocity"
  vvel_azi@units        = "m s~S~-1~N~"

; Create time array for later plotting
  time_int = new(numTIMES,"integer")
  time_int = ispan(0,119,1) ; ispan(ts,tf,1)
  time_int!0 = "time"

;===============================
; Start loop over input files
;===============================

  c  = addfile(fili_c,"r")             ; Read in 'pc' stream [u,v,vort]
  d  = addfile(fili_d,"r")             ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

  ct = 0			       ; Counter variable (time)

;==================================
; Get the variables we will need
;==================================

  do it = ts, tf                   ; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it+1,0)
    t1    = lat_arr(cn,it+1,1)
    n0    = lon_arr(cn,it+1,0)
    n1    = lon_arr(cn,it+1,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Read in basic variables from 'pc' stream
    lon1C  = c->longitude_1({n0:n1})              ; longitude ['d1' grid points]
    lat1C  = c->latitude_1({t0:t1})               ; latitude  ['d0' grid points]
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:)                      ; 63 model ('theta') levels

    u_varname = "x-wind"
    v_varname = "y-wind"

 ; Read in horizontal wind components
 ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
   if (lay .eq. 1) then
    u0  = c->$u_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})     ; Zonal wind (layer)
    v0  = c->$v_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})     ; Meridional wind (layer)
   else
    u0  = c->$u_varname$(it,mlev0,{t0:t1},{n0:n1})           ; Zonal wind (single level)
    v0  = c->$v_varname$(it,mlev0,{t0:t1},{n0:n1})           ; Meridional wind (single level)
   end if

;===================================
; Calculate storm-relative winds
;===================================

    u = u0 - u_cyc(cn,it-1)
    v = v0 - v_cyc(cn,it-1)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

;==================================================================
; Calculate relative vorticity using centered finite differences
;==================================================================

  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u, vort)
    vort@units = "s~S~-1~N~"
    vort@name = "Relative vorticity on model levels"

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(latC(0:dsize(0)-1))

    if (lay .eq. 0) then
     f0   = conform_dims(dimsizes(vort),f1,1)
    elseif (lay .eq. 1) then
     f0   = conform_dims(dimsizes(vort),f1,2)
    end if

    copy_VarMeta(vort,f0)
    copy_VarCoords(vort,f0)
    f0@description = "Coriolis parameter"
    f0@name        = "Coriolis parameter"

    avo  = vort + f0
    copy_VarMeta(vort,avo)
    copy_VarAtts(vort,avo)
    copy_VarCoords(vort,avo)
    avo@description = "Absolute vorticity"
    avo@name        = "Absolute vorticity on model levels"
    avo@long_name   = "Absolute vorticity"

;=====================================================
; Read in variables from 'pd' stream (model levels)
;=====================================================

  ; Read in basic variables from 'pd' stream
    lon1D  = d->longitude_1({n0:n1})                   ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                    ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                     ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                      ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)                           ; 63 model ('rho') levels

    hy     = hybD / 1000 ; Model level values in km
    hy0    = sprintf("%0.0f",hybD(mlev0))
    print("Working on model level: "+hy0+ " m")

    if (lay .eq. 1) then
     hy1    = sprintf("%0.0f",hybD(mlev1))
    end if

  ; Read in temperature, pressure, theta, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]

    if (lay .eq. 1) then
     t      = d->temp(it-2,mlev0:mlev1,{t0:t1},{n0:n1})        ; Temperature (K)
     p      = d->p(it-2,mlev0:mlev1,{t0:t1},{n0:n1})           ; Pressure (Pa)
     th     = d->theta(it-2,mlev0:mlev1,{t0:t1},{n0:n1})       ; Potential temperature (K)
     w      = d->dz_dt(it-2,mlev0:mlev1,{t0:t1},{n0:n1})       ; Vertical velocity (m s-1)
     z      = d->ht(it-2,mlev0:mlev1,{t0:t1},{n0:n1})          ; Geopotential height (m)
     pv	   = d->field83(it-2,mlev0:mlev1,{t0:t1},{n0:n1})      ; Potential vorticity (PVU)
    elseif (lay .eq. 0) then
     t      = d->temp(it-2,mlev0,{t0:t1},{n0:n1})	       ; Temperature (K)
     p      = d->p(it-2,mlev0,{t0:t1},{n0:n1})                 ; Pressure (Pa)
     th     = d->theta(it-2,mlev0,{t0:t1},{n0:n1})             ; Potential temperature (K)
     w      = d->dz_dt(it-2,mlev0,{t0:t1},{n0:n1})             ; Vertical velocity (m s-1)
     z      = d->ht(it-2,mlev0,{t0:t1},{n0:n1})                ; Geopotential height (m)
     pv     = d->field83(it-2,mlev0,{t0:t1},{n0:n1})           ; Potential vorticity (PVU)
    end if

  ; Convert units of pressure and geopotential height before continuing
    p  = p / 100           ; Convert to 'hPa'
    z  = z / 10		   ; Convert to 'dam'
    pv = pv * (10 ^ 6)	   ; Convert to 'PVU' [10 ^ 6 K m2 s-1 kg-1]

  ; Add metadata
    p@units  = "hPa"
    z@units  = "dam"
    pv@units = "PVU"

    dx0      = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0      = latD(1) - latD(0) ; Grid spacing (latitude)

;=======================================
; Calculate storm centre for 'calc=2'
;=======================================

  ; Model level closest to 1 km AGL (980 m)
    ilev = 1

  ; Read in relative vorticity on chosen model level
    vort_plane0 = vort(ilev,:,:)

  ; Calculate location of vorticity maximum (storm centre) on chosen model level
    vort_max0               = max(vort_plane0)
    dims0                   = dimsizes(vort_plane0)
    vort1d0                 = ndtooned(vort_plane0)
    inds0                   = ind_resolve(maxind(vort1d0),dims0)
    delete([/vort_plane0, vort_max0, dims0, vort1d0/])

;===========================================
; Analyse variables on chosen model level
;===========================================

    if (lay .eq. 1) then
     u_plane       = dim_avg_n_Wrap(u(:,:,:),0)
     v_plane       = dim_avg_n_Wrap(v(:,:,:),0)
     vort_plane    = dim_avg_n_Wrap(vort(:,:,:),0)
     prs_plane     = dim_avg_n_Wrap(p(:,:,:),0)
     pv_plane	   = dim_avg_n_Wrap(pv(:,:,:),0)
     w_plane	   = dim_avg_n_Wrap(w(:,:,:),0)
    elseif (lay .eq. 0) then
     u_plane       = u(:,:)
     v_plane       = v(:,:)
     vort_plane    = vort(:,:)
     prs_plane     = p(:,:)
     pv_plane      = pv(:,:)
     w_plane       = w(:,:)
    end if

  ; Max and min pressure on height level
    pmin           = toint(min(prs_plane))
    pmax           = toint(max(prs_plane))

  ; Calculate location of vorticity maximum (storm centre) on model level
    vort_max       = max(vort_plane)
    vort_smth      = smth9_Wrap(vort_plane, 0.5, 0.5, True)
    vort_max_smth  = max(vort_smth)

  ; Find index of relative vorticity maximum
    dims           = dimsizes(vort_plane)
    vort1d         = ndtooned(vort_plane)
    inds           = ind_resolve(maxind(vort1d),dims)
    vort1d_sm      = ndtooned(vort_smth)
    inds_sm        = ind_resolve(maxind(vort1d_sm),dims)

   ;==================================================================
   ; Either calculate the centre, or read in pre-calculated values
   ;==================================================================

   if (calc .eq. 0) then			; Find centre on each model level

   ; Find latitude and longitude of storm centre using vorticity calculation above
     lat_max 	   = t0 + (dy0 * inds(0,0))
     lon_max 	   = n0 + (dx0 * inds(0,1))

     print("ct = "+ct+" ; it = "+it)
     print("lat_max = "+centre(cn,it-1,0)+" ; lon_max = "+centre(cn,it-1,1))

   ; Difference between 'lat_max'/'lon_max' and values read in using NCL function
   ; 'it+1' results from difference in number of times in 'pb' and 'pc' stream arrays
     dlat           = abs(lat_max - centre(cn,it-1,0))
     dlon           = abs(lon_max - centre(cn,it-1,1))

   ; If our calculations of latitude and longitude deviate too much from the values
   ; we read in using the NCL function ('st_rm'), then calculate a new centre
   ; using the minimum pressure. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed.

     if (dlat .gt. 0.75 .or. dlon .gt. 0.75) then

    ; 'it-1' because 'pc' stream data is     our base, and 'centre' has 1 element fewer
      lt1 = centre(cn,it-1,0) - 0.5
      lt2 = centre(cn,it-1,0) + 0.5
      ln1 = centre(cn,it-1,1) - 0.5
      ln2 = centre(cn,it-1,1) + 0.5

      prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid
      dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid
      prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array
      inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum

      lat_max    = lt1 + (dy0 * inds_p0(0,0))
      lon_max    = ln1 + (dx0 * inds_p0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

     end if

   elseif (calc .eq. 1) then			; Find centre using offline script

     lat_max     = centre(cn,it-1,0)
     lon_max     = centre(cn,it-1,1)

   elseif (calc .eq. 2) then			; Find centre using 980 m vorticity maximum

     lat_max     = t0 + (dy0 * inds0(0,0))
     lon_max     = n0 + (dx0 * inds0(0,1))
     print("Centre (980 m): "+lat_max+" degrees N, "+lon_max+" degrees E")

   end if

   ;==============================================================
   ; Create arrays for later switch to cylindrical coordinates
   ;==============================================================

   ;===============
   ; 'pc' stream
   ;===============

   ; Indices of lat/lon grid for 'u' and 'v' in 'pc' stream
     latC_plane  = ind(t0.le.latC.and.latC.le.t1)
     lonC_plane  = ind(n0.le.lonC.and.lonC.le.n1)

   ; Array sizes
     latC_size   = dimsizes(latC_plane)
     lonC_size   = dimsizes(lonC_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     latC_sub    = latC(latC_plane(0):latC_plane(latC_size-1))
     lonC_sub    = lonC(lonC_plane(0):lonC_plane(lonC_size-1))

   ;===============
   ; 'pd' stream
   ;===============

   ; Indices of lat/lon grid for 'str_xz' in 'pd' stream
     latD_plane  = ind(t0.le.latD.and.latD.le.t1)
     lonD_plane  = ind(n0.le.lonD.and.lonD.le.n1)

   ; Array sizes
     latD_size   = dimsizes(latD_plane)
     lonD_size   = dimsizes(lonD_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     latD_sub    = latD(latD_plane(0):latD_plane(latD_size-1))
     lonD_sub    = lonD(lonD_plane(0):lonD_plane(lonD_size-1))

   ; Indices of lat/lon grid for 'str_yz' in 'pd' stream
     lat1D_plane = ind(t0.le.lat1D.and.lat1D.le.t1)
     lon1D_plane = ind(n0.le.lon1D.and.lon1D.le.n1)

   ; Array sizes
     lat1D_size  = dimsizes(lat1D_plane)
     lon1D_size  = dimsizes(lon1D_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat1D_sub   = lat1D(lat1D_plane(0):lat1D_plane(lat1D_size-1))
     lon1D_sub   = lon1D(lon1D_plane(0):lon1D_plane(lon1D_size-1))

     delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

;===========================================
; Switch to cylindrical coordinate system
;===========================================

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing variables we want to plot
     u_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal wind
     v_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional wind
     spd_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vector wind (for RMW)
     pv_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential vorticity
     w_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
     v_rad   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial wind
     v_tan   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential wind

   ; Add variable metadata
     copy_VarAtts(u_plane,u_int)
     u_int@description = "Zonal wind"
     u_int@units       = "m s~S~-1~N~"

     copy_VarAtts(v_plane,v_int)
     v_int@description = "Meridional wind"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(w_plane,w_int)
     v_int@description = "Vertical velocity"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(pv_plane,pv_int)
     pv_int@description = "Potential vorticity"
     pv_int@units       = "10~S~-6~N~K kg~S~-1~N~ m~S~2~N~ s~S~-1~N~"

     v_rad@description = "Radial wind"
     v_rad@units       = "m s~S~-1~N~"

     v_tan@description = "Tangential wind"
     v_tan@units       = "m s~S~-1~N~"

     xpos!0    = "rad"
     xpos!1    = "azi"
     xpos&rad  = radii
     xpos&azi  = thetas

     ypos!0    = "rad"
     ypos!1    = "azi"
     ypos&rad  = radii
     ypos&azi  = thetas

     u_int!0   = "rad"
     u_int!1   = "azi"
     u_int&rad = radii
     u_int&azi = thetas

     v_int!0   = "rad"
     v_int!1   = "azi"
     v_int&rad = radii
     v_int&azi = thetas

     spd_int!0   = "rad"
     spd_int!1   = "azi"
     spd_int&rad = radii
     spd_int&azi = thetas

     w_int!0   = "rad"
     w_int!1   = "azi"
     w_int&rad = radii
     w_int&azi = thetas

     pv_int!0   = "rad"
     pv_int!1   = "azi"
     pv_int&rad = radii
     pv_int&azi = thetas

     v_tan!0   = "rad"
     v_tan!1   = "azi"
     v_tan&rad = radii
     v_tan&azi = thetas

     v_rad!0   = "rad"
     v_rad!1   = "azi"
     v_rad&rad = radii
     v_rad&azi = thetas

   ; Latitude/longitude (small grid)

     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
         theta = thetas_rad(iang)
           xpos(irad,iang) = lon_max + r*cos(theta)
           ypos(irad,iang) = lat_max + r*sin(theta)
       end do
     end do

   ; Fill new arrays with interpolated values on the cylindrical grid
   ; Which latitude/longitude grids are each of these variables defined on?
   ; u_plane     -- latC,  lonC
   ; v_plane     -- latC,  lonC

   ; w_plane     -- latD,  lon1D
   ; pv_plane	 -- latD,  lon1D

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lonC_sub,latC_sub,u_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lonC_sub,latC_sub,v_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1D_sub,latD_sub,w_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      pv_int(:,iang)   = linint2_points(lon1D_sub,latD_sub,pv_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
    end do

  ; Calculate vector windspeed (from 'u' and 'v')
    spd_int = sqrt((u_int ^ 2) + (v_int ^ 2))
    copy_VarCoords(u_int,spd_int)
    copy_VarMeta(u_int,spd_int)
    spd_int@units = "m s~S~-1~N~"

  ; Calculate radial and tangential windspeed at each point
    do iang = 0, dimsizes(thetas)-1
      theta = thetas_rad(iang)

	ra = u_int(:,iang)*cos(theta)
	rb = v_int(:,iang)*sin(theta)

	ta = -u_int(:,iang)*sin(theta)
	tb = v_int(:,iang)*cos(theta)

        v_rad(:,iang) = ra + rb
        v_tan(:,iang) = ta + tb

    end do

  ;===================================================
  ; Average variables azimuthally around the storm
  ;===================================================

    do irad = 0, dimsizes(radii)-1
      r = radii(irad)
      vtan_azi(irad) = avg(v_tan(irad,:))        ; Tangential wind
      vrad_azi(irad) = avg(v_rad(irad,:))	 ; Radial wind
      spd_azi(irad)  = avg(spd_int(irad,:))	 ; Vector wind (for RMW)
      pvrt_azi(irad) = avg(pv_int(irad,:))	 ; Potential vorticity
      vvel_azi(irad) = avg(w_int(irad,:))	 ; Vertical velocity
    end do

  ; Input azimuthal averages into larger arrays before plotting
    vtan_plot0(:,it-1)  = vtan_azi(:)            ; Tangential wind
    vrad_plot0(:,it-1)  = vrad_azi(:)	    	 ; Radial wind
    pvrt_plot0(:,it-1)  = pvrt_azi(:)	    	 ; PV
    vvel_plot0(:,it-1)  = vvel_azi(:)	    	 ; Vertical velocity

  ; Tidy up
    delete([/u_int,v_int,w_int,pv_int,v_tan,v_rad/])
    delete([/u_plane,v_plane,w_plane,vort_plane,prs_plane,pv_plane/])

    ct  = ct + 1		    ; Counter variable (time)

  ; Tidy up (avoid dimension size errors)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort/])
    delete([/lonD,latD,lon1D,lat1D,t,p,w,z,pv/])

   end do     ; End time loop (do it = ts, tf)

;=====================================================
; Calculate radial and tangential wind tendencies
;=====================================================

; Initialise new arrays
  vtan_tend0 	  = vtan_plot0
  vrad_tend0 	  = vrad_plot0

  vtan_tend0!0    = "rad"
  vtan_tend0&rad  = radii
  vtan_tend0!1    = "time"
  vtan_tend0&time = time_int

  vrad_tend0!0    = "rad"
  vrad_tend0&rad  = radii
  vrad_tend0!1    = "time"
  vrad_tend0&time = time_int

  print("Calculating radial and tangential wind tendencies...")

  do it = ts, tf

   tp1  = min((/it+1,tf/))
   tm1  = max((/it-1,ts/))

   vtan_tend0(:,it) = (vtan_plot0(:,tp1) - vtan_plot0(:,tm1)) / 2 ; Units [m s-1 h-1]
   vrad_tend0(:,it) = (vrad_plot0(:,tp1) - vrad_plot0(:,tm1)) / 2 ; Units [m s-1 h-1]

  end do     ; End time loop (do it = ts, tf)

;======================================
; Reorder dimensions before plotting
;======================================

 ; Arrays now ordered [TIME * RADIUS] rather than [RADIUS * TIME]
   vtan_tend = vtan_tend0(time|:,rad|:)
   vrad_tend = vrad_tend0(time|:,rad|:)

   vtan_plot = vtan_plot0(time|:,rad|:)
   vrad_plot = vrad_plot0(time|:,rad|:)
   vvel_plot = vvel_plot0(time|:,rad|:)
   pvrt_plot = pvrt_plot0(time|:,rad|:)

;=====================================================
; Create a panel plot of Figs. 6a and 6b from N11
;=====================================================

  ; Helen Edit file location
  ; Output file location and type
    if (lay .eq. 1) then
     output = "/nfs/see-fs-01_users/earhbu/TropicalCyclones/NCL_out/cyl_coords/hov_uv_"+dat+"_"+ens0+"_"+hy0+"_"+hy1+\
     	      "_"+cn0
    elseif (lay .eq. 0) then
     output = "/nfs/see-fs-01_users/earhbu/TropicalCyclones/NCL_out/cyl_coords/hov_uv_"+dat+"_"+ens0+"_"+hy0+"_"+cn0
    end if
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot and define colourmap
    panel = new(3,graphic)
    gsn_define_colormap(wks,"n11")

;==========================
; Options for plotting
;==========================

  ; Tangential velocity
    opts_vtan                              = True
    opts_vtan@cnFillOn                     = True
    opts_vtan@cnLineLabelInterval          = 2.0
    opts_vtan@cnLineLabelFontHeightF       = 0.012
    opts_vtan@cnLineLabelBackgroundColor   = "transparent"
    opts_vtan@cnLineLabelPlacementMode     = "constant"
    opts_vtan@cnLinesOn                    = False ; Contour lines off
    opts_vtan@cnInfoLabelOn                = False
    opts_vtan@cnLevelSelectionMode         = "ExplicitLevels"

    if (tclr .eq. 0) then
     opts_vtan@cnLevels                   = (/1.0, 2.0, 2.5, 3.0, 3.5, 4.0, \
                                              4.5, 5.0, 6.0, 7.0/)
     opts_vtan@cnFillColors               = (/0,2,3,4,5,6,7,8,9,10,11/)
    elseif (tclr .eq. 1) then
     opts_vtan@cnLevels                   = (/0.5, 1.5, 2.0, 2.5, 3.0, 3.5, \
                                              4.0, 4.5, 5.0, 6.0/)
     opts_vtan@cnFillColors               = (/0,2,3,4,5,6,7,8,9,10,11/)
    elseif (tclr .eq. 2) then
     opts_vtan@cnLevels                   = (/0.5, 1.0, 1.5, 2.0, 2.5, 3.0, \
                                              3.5, 4.0, 4.5, 5.0/)
     opts_vtan@cnFillColors               = (/0,2,3,4,5,6,7,8,9,10,11/)
    end if

    opts_vtan@gsnPaperOrientation          = "landscape"
    opts_vtan@tiMainString                 = ""
    opts_vtan@tiMainFontHeightF            = 0.0125
    opts_vtan@gsnLeftString                = ""
    opts_vtan@gsnRightString               = ""
    opts_vtan@gsnMaximize                  = True
    opts_vtan@lbLabelBarOn                 = True                  ; Labelbar on/off
    opts_vtan@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
    opts_vtan@pmLabelBarWidthF             = 0.32                  ; Labelbar width
    opts_vtan@pmLabelBarHeightF            = 0.08                  ; Labelbar height
    opts_vtan@lbLabelFontHeightF           = 0.0125                ; Labelbar font size
    opts_vtan@lbLabelFont                  = "Helvetica"           ; Labelbar font
    opts_vtan@lbPerimOn                    = False                 ; Perimeter on/off
    opts_vtan@gsnDraw                      = False                 ; Do not draw plot
    opts_vtan@gsnFrame                     = False                 ; Do not advance frame
    opts_vtan@vpWidthF                     = 0.30                  ; Match width in N11
    opts_vtan@vpHeightF                    = 0.50                  ; Match height in N11

  ; Additional plotting resources
    opts_vtan@tiYAxisString                = "Time (h)"
    opts_vtan@tiXAxisString                = "Radius (km)"
    opts_vtan@tiXAxisFontHeightF           = 0.018
    opts_vtan@trYLog                       = False
    opts_vtan@trYMinF                      = ts      ; Starting point (T+...) for plot
    opts_vtan@trYMaxF                      = tf-2    ; Ending point (T+...) for plot
    opts_vtan@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_vtan@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_vtan@gsnMaximize                  = True    ; Maximise plot size
    opts_vtan@gsnAddCyclic                 = False
    opts_vtan@tmXBMode                     = "Explicit"                ; Choose tm mode
    opts_vtan@tmXBValues                   = (/0.0,0.5,1.0,1.5,2.0/)   ; Tm positions
    opts_vtan@tmXBLabels                   = (/"0","50","100",\
    					       "150","200"/)	       ; Tm values
    opts_vtan@tmXBMinorOn                  = False

    opts_vtan@pmLabelBarOrthogonalPosF     = 0.10    ; Move lb up/down (higher/lower)
    opts_vtan@tiXAxisOffsetYF              = 0.10    ; Move x-axis title up/down (h/l)

    opts_vtan@trXMinF                      = 0.0                       ; Set x-axis min
    opts_vtan@trXMaxF                      = ar0                       ; Set x-axis max
    opts_vtan@tmYLMode                     = "Manual"
    opts_vtan@tmYLTickStartF               = opts_vtan@trYMinF         ; See L730
    opts_vtan@tmYLTickEndF                 = opts_vtan@trYMaxF         ; See L731
    opts_vtan@tmYLTickSpacingF             = 6
    opts_vtan@tmYLMinorOn                  = False

  ; Resources to overlay contours of mean tangential wind
    opts_cont				   = True

    opts_cont@vpWidthF                     = 0.30                  ; Match width in N11
    opts_cont@vpHeightF                    = 0.50                  ; Match height in N11

    opts_cont@cnFillOn                     = False       ; Filled contours on/off
    opts_cont@cnLineColor                  = "blue3"     ; Contour line colour
    opts_cont@cnLevelSelectionMode         = "ExplicitLevels"
    opts_cont@cnLevels                     = (/10,20,30,40,50,60,70/)

    opts_cont@tiYAxisString                = "Time (h)"
    opts_cont@tiXAxisString                = "Radius (km)"
    opts_cont@tiXAxisFontHeightF           = 0.018
    opts_cont@trYLog                       = False
    opts_cont@trYMinF                      = ts          ; Starting point (T+...) for plot
    opts_cont@trYMaxF                      = tf-2        ; Ending point (T+...) for plot
    opts_cont@tmXTOn                       = "False" 	 ; Turn off top x-axis TM
    opts_cont@tmYROn                       = "False" 	 ; Turn off right y-axis TM
    opts_cont@gsnMaximize                  = True    	 ; Maximise plot size
    opts_cont@gsnAddCyclic                 = False
    opts_cont@tmXBMode                     = "Explicit"                ; Choose tm mode
    opts_cont@tmXBValues                   = (/0.0,0.5,1.0,1.5/)       ; Tm positions
    opts_cont@tmXBLabels                   = (/"0","50","100","150"/)  ; Tm values
    opts_cont@tmXBMinorOn                  = False

    opts_cont@pmLabelBarOrthogonalPosF     = 0.10    ; Move lb up/down (higher/lower)
    opts_cont@tiXAxisOffsetYF              = 0.10    ; Move x-axis title up/down (h/l)

    opts_cont@trXMinF                      = 0.0                       ; Set x-axis min
    opts_cont@trXMaxF                      = 1.0                       ; Set x-axis max

    opts_cont@cnLineThicknessF             = 2.0         ; Contour line thickness
    opts_cont@cnInfoLabelOn                = False       ; Info label on/off
    opts_cont@gsnDraw                      = False       ; Do not draw the plot
    opts_cont@gsnFrame                     = False       ; Do not advance the frame
    opts_cont@cnLineLabelsOn               = True        ; Turn line labels on/off
    opts_cont@cnLineLabelFontColor         = "blue3"     ; Line label colour
    opts_cont@cnLineLabelFontHeightF       = 0.01        ; Default 0.006 for vpWidth=0.30
    opts_cont@cnLineLabelFontThicknessF    = 2.0         ; Line label thickness
    opts_cont@tiMainString                 = ""
    opts_cont@gsnLeftString                = ""
    opts_cont@gsnRightString               = ""

  ; Resources to overlay zero line
    opts_z                             = True
    opts_z@cnFillOn                    = False
    opts_z@cnLineColor                 = "gray82"
    opts_z@cnLevels                    = (/-1.0, 0.0, 1.0/)
    opts_z@cnLineThicknessF            = 2.0
    opts_z@cnInfoLabelOn               = False
    opts_z@gsnDraw                     = False       ; Do not draw the plot
    opts_z@gsnFrame                    = False       ; Do no advance the frame
    opts_z@cnLineLabelsOn              = False       ; Turn off line labels
    opts_z@tiMainString                = ""
    opts_z@gsnLeftString               = ""
    opts_z@gsnRightString              = ""

  ; Radial velocity resources
    opts_vrad                          = opts_vtan
    delete([/opts_vrad@cnLevels,opts_vrad@cnFillColors/])

    if (rclr .eq. 0) then
     opts_vrad@cnLevels                     = (/0.0, -0.5, -1.0, -1.5, -2.0, \
                                                -2.5, -3.0, -4.0, -5.0/)
     opts_vrad@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    elseif (rclr .eq. 1) then
     opts_vrad@cnLevels                     = (/0.0, -0.25, -0.5, -1.0, -1.5, \
                                                -2.0, -2.5, -3.0, -4.0/)
     opts_vrad@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    elseif (rclr .eq. 2) then
     opts_vrad@cnLevels                     = (/0.0, -0.25, -0.5, -0.75, -1.0, \
                                                -1.5, -2.0, -2.5, -3.0/)
     opts_vrad@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    end if

  ; Vertical velocity resources
    opts_vvel                                = opts_vtan
    delete([/opts_vvel@cnLevels,opts_vvel@cnFillColors/])

    if (wclr .eq. 0) then
     opts_vvel@cnLevels                      = (/0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                                 0.7, 0.8, 0.9, 1.0/)
     opts_vvel@cnFillColors                  = (/0,2,3,4,5,6,7,8,9,10,11,12/)
    elseif (wclr .eq. 1) then
     opts_vvel@cnLevels                      = (/0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, \
                                                 0.50, 0.60, 0.70, 0.80/)
     opts_vvel@cnFillColors                  = (/0,2,3,4,5,6,7,8,9,10,11,12/)

    elseif (wclr .eq. 2) then
     opts_vvel@cnLevels                      = (/0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, \
                                                 0.40, 0.45, 0.50, 0.60/)
     opts_vvel@cnFillColors                  = (/0,2,3,4,5,6,7,8,9,10,11,12/)
    end if

  ; Radius of maximum wind (RMW) resources
    opts_rmw                               = True
    opts_rmw@xyLineColors                  = "red3"
    opts_rmw@xyLineThicknesses             = 3.0
    opts_rmw@xyDashPatterns                = 0.0
    opts_rmw@gsnDraw                       = False       ; Do not draw the plot
    opts_rmw@gsnFrame                      = False       ; Do not advance the frame
    opts_rmw@tiMainString                  = ""
    opts_rmw@gsnLeftString                 = ""
    opts_rmw@gsnRightString                = ""

;========================
; Construct panel plot
;========================

  ; Read in RMW data from text file (first 3 values equal to zero)
    if (rmw .eq. 1) then
     spd_plot = asciiread("/nfs/a37/earshar/ncl_scripts/nepartak/archer_march2018/text/mlev_"+dat+\
                          "_"+ens0+"_rmw_calc1.txt",(/numTIMES/),"float")
    elseif (rmw .eq. 0) then
     spd_plot = asciiread("/nfs/a37/earshar/ncl_scripts/nepartak/archer_march2018/text/"+dat+\
                          "_"+ens0+"_rmw_calc1.txt",(/numTIMES/),"float")
    end if

  ; Tangential wind (LHS panel)
    vtan0 = gsn_csm_contour(wks,vtan_tend,opts_vtan)    ; Tangential wind tendency (filled)
    zero0 = gsn_csm_contour(wks,vtan_tend,opts_z)    	; Zero line (vtan tendency)
    zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
    cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)    ; Tangential wind
    rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
    overlay(vtan0,zero)                                 ; Overlay vtan tendency zero line
    overlay(vtan0,rmw0)				        ; Overlay RMW
    overlay(vtan0,cont0)                                ; Overlay tangential wind
    panel(0) = vtan0                                    ; Plot first panel
    delete([/zero0,zero/])                              ; Tidy up

  ; Radial wind (middle panel)
    vrad0 = gsn_csm_contour(wks,vrad_tend,opts_vrad)    ; Radial wind tendency
    cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)    ; Tangential wind
    rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
    overlay(vrad0,cont0)                                ; Overlay tangential wind
    overlay(vrad0,rmw0)                                 ; Overlay RMW
    panel(1) = vrad0                                    ; Plot second panel

  ; Vertical velocity (RHS panel)
    vvel0 = gsn_csm_contour(wks,vvel_plot,opts_vvel)    ; Vertical velocity
    zero0 = gsn_csm_contour(wks,vvel_plot,opts_z)       ; Zero line (vertical velocity)
    zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
    rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
    overlay(vvel0,zero)                                 ; Overlay vvel zero line
    overlay(vvel0,rmw0)                                 ; Overlay RMW
    panel(2) = vvel0                                    ; Plot third panel
    delete([/zero0,zero/])                              ; Tidy up

;=============================================
; Finally, put everything together and plot
;=============================================

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn on panel labelbar
    optsP@txString              = "Simulation "+ens0+" ("+dat+")"

    optsP@gsnPanelFigureStrings = (/"a) vtan_tend","b) vrad_tend","c) vvel"/)

    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Reduce label size (default 0.01)

    gsn_panel(wks,panel,(/1,3/),optsP)             ; Draw as a single plot
    frame(wks)

end
