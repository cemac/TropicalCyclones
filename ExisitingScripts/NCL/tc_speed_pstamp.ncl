; Script to produce postage stamp plot of horizontal windspeed for ensemble members

; Run using

; ncl dat=\"03T12\" dist=3.0 plev=3 opt=\"pdf\" cn0=\"geo_sm\" calc=1
; ts0=24 tf0=30 tc_speed_pstamp_cemac.ncl

; 'int'  = interval in the time loop
; 'dist' = size of box following storm (degrees)
; 'opt'  = output file format ("pdf","x11",etc)
; 'plev' = pressure level index (2,5,7,etc)
; 'dat'  = initialisation time string ("02T12", "03T00", ...)
; 'cn0'  = method of storm tracking ["slp", "geo", "geo_sm", "vort"] - default is "geo_sm"
; 'ts0'  = starting time index (3-119)
; 'tf0'  = ending time index (3-119)

; Load main NCL libraries ($NCARG_ROOT = /apps/applications/ncl/6.5.0/1/default)
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/nfs/a37/earshar/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

; Also load additional user-written functions
load "ncl_func/st_rm.ncl"
load "ncl_func/nc_times.ncl"
load "ncl_func/nc_levs.ncl"

begin

; Exit script if missing times are chosen
  if (ts0 .lt. 3) then
   print("Starting time index must be 3 or greater (data missing)! Please try again...")
   exit()
  end if

;=============================================
; Define set of simulations we want to plot
;=============================================

  ens_arr     = (/"em00","em01","em02","em03","em04","em05",\
                  "em06","em07","em08","em09","em10","em11"/)

  numSIMS     = dimsizes(ens_arr)

; Get file information before starting loop below
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens_arr(0)                 ; Directory path
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens_arr(0)        ; Main file name
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")   ; File type
  fili_b      = diri+"/"+fili_b1

; Call user-defined function and get number of times in file
  storm_rel   = st_rm(dat, ens_arr(0), "$ar/text/", dist, 0)
  centre      = storm_rel[6]             ; Storm track (position) array
  cdims	      = dimsizes(centre)         ; [4 * 120 * 2]

; Create arrays to hold values output from 'st_rm' multiple times (12 panels)
  centre_all  = new((/numSIMS,cdims(0),cdims(1),cdims(2)/),float)
  lat_all     = new((/numSIMS,cdims(0),cdims(1),cdims(2)/),float)
  lon_all     = new((/numSIMS,cdims(0),cdims(1),cdims(2)/),float)
  uc_all      = new((/numSIMS,cdims(0),cdims(1)/),float)
  vc_all      = new((/numSIMS,cdims(0),cdims(1)/),float)

;=========================================================================
; Also create array to hold windspeed data for all simulations in panel
;=========================================================================

; Create arrays for data and output information
  llbox                = toint(dist*50)           ; Calculate domain size (grid points)
  dsize                = (/llbox,llbox/)          ; Domain size (grid points)
  d0                   = dsize(0)                 ; Size of dimension 0
  d1                   = dsize(1)                 ; Size of dimension 1

; Create 4D array --> [12 * 120 * 150 * 150]
  speed_all   	       = new((/numSIMS,cdims(1),d0,d1/),float)

; Also create 3D arrays to hold values of latitude and longitude for each simulation, at each time
  pos_lat              = new((/numSIMS,cdims(1),d0/),float)
  pos_lon              = new((/numSIMS,cdims(1),d1/),float)

; Tidy up before main loop
  delete([/diri,fili_prefix,fili_b1,fili_b,storm_rel,centre,cdims/])

;================================
; Loop over chosen simulations
;================================

  do sm        = 0, numSIMS-1

;===============================================
; Find 'pb' stream data (pressure level data)
;===============================================

   print("Reading in storm track data for "+ens_arr(sm)+" ("+dat+")")

   diri        = "/nfs/a319/earshar/"+dat+"/"+ens_arr(sm)
   fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens_arr(sm)
   fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb.nc")
   fili_b      = diri+"/"+fili_b1

;=============================================================
; Calculate number of times in file using built-in function
;=============================================================

 ; 'pb' stream
   numINFO_b   = nc_times(fili_b)
   numINPUT_b  = numINFO_b[0]            ; Number of times
   time_b      = numINFO_b[1]            ; Array of times
   print("File contains "+numINPUT_b+" times")

;===============================================
; Calculate number of vertical levels in file
;===============================================

 ; 'pb' stream
   nLEVS_b     = nc_levs(fili_b)
   numLEVS_b   = nLEVS_b[1]              ; Number of vertical levels
   levs_b      = nLEVS_b[0]              ; Array of vertical levels

;============================================================
; Calculate storm motion using built-in function ('st_rm')
;============================================================

 ; 'dat'       = initialisation time (02T12, 03T00, ...)
 ; 'ens0'      = ensemble simulation (em00, em01, ...)
 ; 'diri'      = path to input files (see above)
 ; 'dist'      = size of box following storm (degrees)
 ; 'mins'      = analyse 1-h (0) or 5-min (1) data

   storm_rel   	        = st_rm(dat, ens_arr(sm), "$ar/text/", dist, 0)
   numINPUT_a  	        = storm_rel[3]             ; Number of times in file
   uc_all(sm,:,:)       = storm_rel[0]             ; Zonal wind
   vc_all(sm,:,:)       = storm_rel[1]             ; Meridional wind
   lat_all(sm,:,:,:)    = storm_rel[4]             ; Latitude array subset (following storm)
   lon_all(sm,:,:,:)    = storm_rel[5]             ; Longitude array subset (following storm)
   centre_all(sm,:,:,:) = storm_rel[6]             ; Storm track (position) array

 ; Increase memory
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
; Create date/time string arrays for all times in file (use 'pb' stream)
;==========================================================================

   month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                     "Jul","Aug","Sep","Oct","Nov","Dec"/)

   utc_date = cd_calendar(time_b, 0)

   year   = tointeger(utc_date(:,0))
   month  = tointeger(utc_date(:,1))
   day    = tointeger(utc_date(:,2))
   hour   = tointeger(utc_date(:,3))
   minute = tointeger(utc_date(:,4))
   second = utc_date(:,5)

 ; Correct for errors in the code (round up value of hour when minutes = 59)
   do ct = 0, numINPUT_b-1
    if (minute(ct).gt.30) then
      hour(ct) = hour(ct)+1
    end if
   end do
   delete(ct)

 ; Create arrays to hold strings for output plots (117 elements)
   date_str  = new(numINPUT_b,string)
   time_str  = new(numINPUT_b,string)
   time_arr  = new(numINPUT_b,string)
   title_arr = new(numINPUT_b,string)

   date_str  = sprinti("%0.2iUTC ", hour) + \
               sprinti("%0.2i ", day) \
             + month_abbr(month)
   time_str  = sprinti("%0.2i UTC ", hour) + \
               sprinti("%0.2i ", day) \
             + month_abbr(month)
   time_arr  = sprinti("%0.2i", day) + \
               month_abbr(month) + \
               "_" + sprinti("%0.2iZ", hour)

;==================================
; Get the variables we will need
;==================================

  b         = addfile(fili_b,"r")               ; Read in 'pb' stream [u,v,w,vort,t,z,rh]

  do it     = ts0, tf0

   print("Working on time: "+time_str(it)+" (it = "+it+")")
   itt = it + 1
   title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

   time  = b->t(it)           ; Time in file
   times = dimsizes(time)

   t0    = lat_all(sm,cn,it,0)
   t1    = lat_all(sm,cn,it,1)
   n0    = lon_all(sm,cn,it,0)
   n1    = lon_all(sm,cn,it,1)
   print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

 ; Create latitude/longitude array subsets to reduce computation time
   lon1  = b->longitude_1({n0:n1})  ; longitude ['d0' grid points]
   lat1  = b->latitude_1({t0:t1})   ; latitude  ['d1' grid points]
   lon   = b->longitude({n0:n1})    ; longitude ['d0' grid points]
   lat   = b->latitude({t0:t1})     ; latitude  ['d1' grid points]

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800, (6) 750, (7) 700, (8) 650,
 ; (9) 600, (10) 500, (11) 400, (12) 300, (13) 250, (14) 200, (15) 150, (16) 100

   u0    = b->u(it-3,plev,{t0:t1},{n0:n1})     ; Zonal velocity (m/s)
   v0    = b->v(it-3,plev,{t0:t1},{n0:n1})     ; Meridional velocity (m/s)
   vort  = b->rvor(it-3,plev,{t0:t1},{n0:n1})  ; Relative vorticity (/s)

   t     = b->temp(it-3,plev,{t0:t1},{n0:n1})  ; Temperature (K)
   z     = b->ht(it-3,plev,{t0:t1},{n0:n1})    ; Geopotential height (m)
   w     = b->dz_dt(it-3,plev,{t0:t1},{n0:n1}) ; Vertical velocity (m/s)

   pres  = b->p                                ; Pressure (hPa)
   plevs = dimsizes(pres)
   pres@units = "hPa"

 ; Convert to units we want
   z     = z/10
   z@units = "dam" ; Convert to decametres
   vort  = vort * (10 ^ -6) ; Convert to /s

   dx0   = lon(1) - lon(0) ; Grid spacing (longitude)
   dy0   = lat(1) - lat(0) ; Grid spacing (latitude)

;==============================================================
; Make sure all arrays are same size (for later calculations)
;==============================================================

 ; 2-D array (6 * 2)
   size  = (/dimsizes(u0),dimsizes(v0),dimsizes(vort),dimsizes(t),dimsizes(z),dimsizes(w)/)

 ; Reduce size of any dimension larger than the rest (u,v,w,t,etc)
   do sz  = 0, 1
    nind = min( size(:,sz) )
    if ( any (size(:,sz) .gt. nind) ) then
     size(:,sz) = nind
    end if
   end do

 ; Resize arrays based upon results of calculations above
   u2    = u0(0:size(0,0)-1,0:size(0,1)-1)
   v2    = v0(0:size(1,0)-1,0:size(1,1)-1)
   vort2 = vort(0:size(2,0)-1,0:size(2,1)-1)
   t2    = t(0:size(3,0)-1,0:size(3,1)-1)
   z2    = z(0:size(4,0)-1,0:size(4,1)-1)
   w2    = w(0:size(5,0)-1,0:size(5,1)-1)
   lat2  = lat(0:nind-1)
   lon2  = lon(0:nind-1)
   delete([/u0,v0,vort,t,z,w,lat,lon/])

 ; Reset names to avoid further code changes
   u0    = u2
   v0    = v2
   vort  = vort2
   t     = t2
   z     = z2
   w     = w2
   lat   = lat2
   lon   = lon2
   delete([/u2,v2,vort2,t2,z2,w2,lat2,lon2/])

;===================================
; Calculate storm-relative winds
;===================================

   u = u0 - vc_all(sm,cn,it)
   v = v0 - vc_all(sm,cn,it)

 ; Add metadata from original horizontal wind arrays
   copy_VarCoords(u0,u)
   copy_VarAtts(u0,u)

   copy_VarCoords(v0,v)
   copy_VarAtts(v0,v)

;=================================
; Pressure levels ['pb' stream]
;=================================

 ; Rename single-layer variables
   u_plane       = u
   v_plane       = v

 ; Rename all other variables
   vort_plane    = vort
   geo_plane     = z
   t_plane       = t

 ; Calculate minimum and maximum geopotential height on pressure level
   zmin          = toint(min(geo_plane))
   zmax          = toint(max(geo_plane))

 ; Locate relative vorticity maximum (storm centre) on pressure level
   vort_max      = max(vort_plane)
   vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
   vort_max_smth = max(vort_smth)

 ; Find index of relative vorticity maximum
   dims          = dimsizes(vort_plane)
   vort1d        = ndtooned(vort_plane)
   inds          = ind_resolve(maxind(vort1d),dims)
   vort1d_sm     = ndtooned(vort_smth)
   inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

;==================================================================
; Either calculate the centre, or read in pre-calculated values
;==================================================================

   if (calc .eq. 0) then

   ; Find latitude and longitude of storm centre using vorticity calculation above
     lat_max0        = t0 + (dy0 * inds(0,0))
     lon_max0        = n0 + (dx0 * inds(0,1))

     print("Centre (vort): "+lat_max0+" degrees N, "+lon_max0+" degrees E")

   ; Create smaller grid around location identified using vorticity maximum
     lt1 = lat_max0 - 0.15
     lt2 = lat_max0 + 0.15
     ln1 = lon_max0 - 0.15
     ln2 = lon_max0 + 0.15

   ; Calculate minimum geopotential height in region immediately surrounding 'storm centre'
     geo_plane0 = geo_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
     geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
     dims_z0    = dimsizes(geo_plane0)                 ; Size of smaller grid
     geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
     inds_z0    = ind_resolve(minind(geo_1d0),dims_z0) ; Find index of minimum

     lat_max0   = lt1 + (dy0 * inds_z0(0,0))
     lon_max0   = ln1 + (dx0 * inds_z0(0,1))

     print("Centre (vort + geo): "+lat_max0+" degrees N, "+lon_max0+" degrees E")
     delete([/geo_plane0,geo_min0,dims_z0,geo_1d0,inds_z0/])

   ; Distance between centre calculated above using relative vorticity, and that calculated offline
     dlat       = abs(lat_max0 - centre_all(sm,cn,it,0))
     dlon       = abs(lon_max0 - centre_all(sm,cn,it,1))

   ; If our storm centre calculations deviate too much (> 0.75º) from the values
   ; we read in using the NCL function 'st_rm', calculate a new centre using
   ; the minimum geopotential height. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed.

     if (dlat .gt. 1.00 .or. dlon .gt. 1.00) then

      lt1 = centre_all(sm,cn,it,0) - 0.25
      lt2 = centre_all(sm,cn,it,0) + 0.25
      ln1 = centre_all(sm,cn,it,1) - 0.25
      ln2 = centre_all(sm,cn,it,1) + 0.25

      geo_plane0 = geo_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
      dims_z0    = dimsizes(geo_plane0)                 ; Size of smaller grid
      geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
      inds_z0    = ind_resolve(minind(geo_1d0),dims_z0) ; Find index of minimum

      lat_max0   = lt1 + (dy0 * inds_z0(0,0))
      lon_max0   = ln1 + (dx0 * inds_z0(0,1))
      print("New centre: "+lat_max0+" degrees N, "+lon_max0+" degrees E")
      delete([/geo_plane0,geo_min0,dims_z0,geo_1d0,inds_z0/])

     end if

   elseif (calc .eq. 1) then

      lat_max0   = centre_all(sm,cn,it,0)
      lon_max0   = centre_all(sm,cn,it,1)
      print("Storm centre ("+ens_arr(sm)+") = ["+lat_max0+" , "+lon_max0+"]")

   end if

;========================================================
; Calculate vector windspeed and write to larger array
;========================================================

   speed_all(sm,it,:,:) = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )

 ; Add metadata
   copy_VarCoords(u_plane,speed_all(sm,it,:,:))
   copy_VarAtts(u_plane,speed_all(sm,it,:,:))

 ; Add latitude and longitude values to larger arrays 'pos_lat' and 'pos_lon'
   pos_lat(sm,it,:) = u_plane&latitude
   pos_lon(sm,it,:) = u_plane&longitude

 ; Tidy up (avoid dimension size errors)
   delete([/lon,lat,lon1,lat1,u0,v0,u,v,vort,t,z,w/])

  end do  ; End time loop (do it = ts0, tf0)

 end do	; End ensemble member loop (do sm = 0, numSIMS-1)

;===============================================
; Produce panel plots for each time interval
;===============================================

 print("Working on panel plots!")
 do it = ts0, tf0

 ; Create object to hold multiple panels
   panel       = new(numSIMS,graphic)

;=================================
; Output file type and location
;=================================

   ; Helen Edit:
   ;output = "/nfs/a37/earshar/nepartak/images/cart_coords/pstamp_uv"+pres(plev)+\
   ;         "_"+dat+"_"+time_arr(it)
   output = "/nfs/see-fs-01_users/earhbu/TropicalCyclones/NCL_out/cart_coords/pstamp_uv"+pres(plev)+\
            "_"+dat+"_"+time_arr(it)
   wks    = gsn_open_wks(opt,output)

;==============================================================
; Load NIMROD colour table for precip duration data
;==============================================================

 ; Load colour table
   gsn_define_colormap(wks,"radar_new")

 ; Change the colour map around & make first colour transparent
   cmap = read_colormap_file("radar_new")
   cmap(0,3) = 0.0

;==========================
; Options for plotting
;==========================

 ; Horizontal windspeed resources
   opts_wind                              = True
   opts_wind@cnFillOn                     = True              ; Turn contour fill on
   opts_wind@cnLinesOn                    = False             ; Contour lines off
   opts_wind@cnLineLabelsOn               = False             ; Contour line labels off
   opts_wind@cnInfoLabelOn                = False             ; Contour label info off
   opts_wind@cnLevelSelectionMode         = "ExplicitLevels"  ; Explicit level selection
   opts_wind@cnLevels                     = (/5.0, 10.0, 15.0, 20.0,\             ; Contour levels
                                              25.0, 30.0, 35.0, 40.0, \
                                              45.0, 50.0, 55.0, 60.0, 65.0/)
   opts_wind@cnFillColors 		  = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14/)  ; Contour colours
   opts_wind@lbLabelBarOn                 = False             ; Turn off labelbar (for now)
   opts_wind@gsnPaperOrientation          = "landscape"

 ; Title resources (no text)
   opts_wind@tiMainString                 = ""
   opts_wind@gsnLeftString                = ""
   opts_wind@gsnRightString               = ""
   opts_wind@tiMainFontHeightF            = 0.0125            ; Title font size

 ; Additional plotting resources
   opts_wind@tmXBMajorLengthF             = 0.005      ; Change tickmark length
   opts_wind@tmXBLabelFontHeightF         = 0.021
   opts_wind@tmYLLabelFontHeightF         = 0.021
   opts_wind@tmXTOn                       = "False"    ; No tickmarks on top x-axis
   opts_wind@tmYROn                       = "False"    ; No tickmarks on right y-axis
   opts_wind@gsnMaximize                  = True       ; Maximise plot size
   opts_wind@gsnAddCyclic                 = False
   opts_wind@gsnDraw                      = False      ; Do not draw the plot
   opts_wind@gsnFrame                     = False      ; Do not advance the frame

;===================================================
; Customise lat/lon labels depending on plot size
;===================================================

 if (dist .le. 2.0) then

 ; Bottom x-axis
   opts_wind@tmXBMode                     = "Explicit"
   opts_wind@tmXBValues                   = fspan(100.0, 150.0, 51)
   xb_labels0                             = sprintf("%0.1f~S~o~N~E",opts_wind@tmXBValues)
   opts_wind@tmXBLabels                   = xb_labels0

 ; Left y-axis
   opts_wind@tmYLMode                     = "Explicit"
   opts_wind@tmYLValues                   = fspan(0.0, 30.0, 31)
   yl_labels0                             = sprintf("%0.1f~S~o~N~E",opts_wind@tmYLValues)
   opts_wind@tmYLLabels                   = yl_labels0

 else

 ; Bottom x-axis
   opts_wind@tmXBMode                     = "Explicit"
   opts_wind@tmXBValues                   = fspan(100.0, 150.0, 26)
   xb_labels0                             = sprintf("%0.1f~S~o~N~E",opts_wind@tmXBValues)
   opts_wind@tmXBLabels                   = xb_labels0

 ; Left y-axis
   opts_wind@tmYLMode                     = "Explicit"
   opts_wind@tmYLValues                   = fspan(0.0, 30.0, 16)
   yl_labels0                             = sprintf("%0.1f~S~o~N~E",opts_wind@tmYLValues)
   opts_wind@tmYLLabels                   = yl_labels0

 end if

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

 ; Plot each simulation using 'gsn_csm_contour' function
   do ii      = 0, dimsizes(panel)-1

  ; Add in correct latitude & longitude information before plotting
    speed_plane = speed_all(ii,it,:,:)          ; Create temporary copy of 'speed_all(mm,tt,:,:)'
    speed_plane&latitude = pos_lat(ii,it,:)     ; Add latitude grid for chosen simulation
    speed_plane&longitude = pos_lon(ii,it,:)    ; Add longitude grid for chosen simulation

  ; Add updated metadata to existing array
    speed_all(ii,it,:,:) = speed_plane

    panel(ii) = gsn_csm_contour(wks, speed_all(ii,it,:,:), opts_wind )

   end do

 ; Create panel plot with all ensemble members
   optsP                       = True
   optsP@gsnFrame              = False    ; Do not advance the frame

 ; Turn on panel labelbar and customise size of final plot
   optsP@gsnPanelLabelBar      = True
   optsP@amJust                = "TopLeft"
   optsP@gsnPanelTop           = 0.98
   optsP@gsnPanelBottom        = 0.02
   optsP@gsnMaximize           = True

 ; Labelbar resources
   optsP@pmLabelBarWidthF      = 0.6
   optsP@pmLabelBarHeightF     = 0.15
   optsP@lbLabelFontHeightF    = 0.01
   optsP@lbPerimOn             = False
   optsP@lbLabelFont           = "Helvetica"
   optsP@lbTitleString         = pres(plev)+" hPa horizontal windspeed (m s~S~-1~N~)"
   optsP@lbTitleFontHeightF    = 0.015
   optsP@lbTitleDirection      = "Across"
   optsP@lbTitlePosition       = "Bottom"
   optsP@lbBoxEndCapStyle      = "TriangleBothEnds"    ; Labelbar end shape

 ; Leave title blank
   optsP@txString              = ""

 ; Panel figure labels
   optsP@gsnPanelFigureStrings = (/"a) em00", "b) em01", "c) em02", "d) em03", \
   			       	   "e) em04", "f) em05", "g) em06", "h) em07", \
                                   "i) em08", "j) em09", "k) em10", "l) em11"/)

 ; Reduce size of labels (default = 0.01)
   optsP@gsnPanelFigureStringsFontHeightF = 0.0099

 ; Produce panel with 3 rows and 4 columns
   gsn_panel(wks,panel,(/3,4/),optsP)
   frame(wks)

  end do     ; End time loop (do it = ts0, tf0)

end
