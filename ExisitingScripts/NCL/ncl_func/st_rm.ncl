;==================================
; Function: st_rm
; Note: calculate storm motion
;==================================

undef("st_rm")
function st_rm(dat:string,ens0:string,diri:string,dist:numeric,mins:numeric)
local ntimes, ntimes_b, centre, dat, ens0, xc, yc, uc, vc, velc, lt_arr, ln_arr

begin

; Number of times in input files
  if (mins .eq. 0) then
   ntimes   = 120			       ; Number of times
   ntimes_b = 117
  elseif (mins .eq. 1) then
   ntimes   = 996
  end if

; Create array to hold storm track (lat/lon) at each time
  if (mins .eq. 0) then
   centre = new((/4,ntimes,2/),float)
  elseif (mins .eq. 1) then
   centre = new((/2,ntimes,2/),float)
  end if

; Initialise with zeros
  centre = 0

; Read in values from chosen simulation
  if (mins .eq. 0) then
   print(diri+"cp_new_"+dat+"_"+ens0+"_lat.txt")  
   centre(0,:,0)   = asciiread(diri+"cp_new_"+dat+"_"+ens0+"_lat.txt",(/ntimes,1/),"float")
   centre(0,:,1)   = asciiread(diri+"cp_new_"+dat+"_"+ens0+"_lon.txt",(/ntimes,1/),"float")
   centre(1,3::,0) = asciiread(diri+"cp_new_"+dat+"_"+ens0+"_lat_rv.txt",(/ntimes_b,1/),"float")
   centre(1,3::,1) = asciiread(diri+"cp_new_"+dat+"_"+ens0+"_lon_rv.txt",(/ntimes_b,1/),"float")
   centre(2,3::,0) = asciiread(diri+"cp_new_"+dat+"_"+ens0+"_lat_z0.txt",(/ntimes_b,1/),"float")
   centre(2,3::,1) = asciiread(diri+"cp_new_"+dat+"_"+ens0+"_lon_z0.txt",(/ntimes_b,1/),"float")
   centre(3,3::,0) = asciiread(diri+"cp_new_"+dat+"_"+ens0+"_lat_z1.txt",(/ntimes_b,1/),"float")
   centre(3,3::,1) = asciiread(diri+"cp_new_"+dat+"_"+ens0+"_lon_z1.txt",(/ntimes_b,1/),"float")
  elseif (mins .eq. 1) then
   centre(0,:,0)   = asciiread(diri+"cp_5min_"+dat+"_"+ens0+"_lat.txt",(/ntimes,1/),"float")
   centre(0,:,1)   = asciiread(diri+"cp_5min_"+dat+"_"+ens0+"_lon.txt",(/ntimes,1/),"float")
   centre(1,:,0)   = asciiread(diri+"cp_5min_"+dat+"_"+ens0+"_latf.txt",(/ntimes,1/),"float")
   centre(1,:,1)   = asciiread(diri+"cp_5min_"+dat+"_"+ens0+"_lonf.txt",(/ntimes,1/),"float")
  end if

  cr = inttochar(10)
  print_clock("Starting storm motion calculation for "+ens0)

; Define arrays to hold values of storm position/motion
  if (mins .eq. 0) then
   xc	  = new((/4,ntimes/),float)
   yc	  = new((/4,ntimes/),float)
   uc	  = new((/4,ntimes/),float)
   vc	  = new((/4,ntimes/),float)
   velc	  = new((/4,ntimes/),float)
  elseif (mins .eq. 1) then
   xc     = new((/2,ntimes/),float)
   yc     = new((/2,ntimes/),float)
   uc     = new((/2,ntimes/),float)
   vc     = new((/2,ntimes/),float)
   velc   = new((/2,ntimes/),float)
  end if

;===============================================================
; Calculate system motion from 1-h cyclone position [lat,lon]
;===============================================================

  if (mins .eq. 0) then

   do ct    = 0, 3
    do iang = 1, ntimes-1
     xc(ct,iang)   = gc_latlon(centre(ct,iang-1,0),centre(ct,iang-1,1),\
                               centre(ct,iang-1,0),centre(ct,iang,1),2,3)
     yc(ct,iang)   = gc_latlon(centre(ct,iang-1,0),centre(ct,iang-1,1),\
                               centre(ct,iang,0),centre(ct,iang-1,1),2,3)

   ; If storm is moving westward, zonal wind is negative
     if ( centre(ct,iang,1) .lt. centre(ct,iang-1,1) ) then
      uc(ct,iang)   = -( xc(ct,iang)/(3600) )
     else
      uc(ct,iang)   = xc(ct,iang)/(3600)
     end if

   ; If storm is moving southward, meridional wind is negative
     if ( centre(ct,iang,0) .lt. centre(ct,iang-1,0) ) then
      vc(ct,iang)   = -( yc(ct,iang)/(3600) )
     else
      vc(ct,iang)   = yc(ct,iang)/(3600)
     end if

   ; Calculate vector wind
     velc(ct,iang) = sqrt( uc(ct,iang) ^ 2.0 + vc(ct,iang) ^ 2.0)

    end do      	; End loop over times
   end do		; End loop over storm tracker methods

 ; Set initial system velocity [0] equal to that at the second time [1]
   xc(:,0)   = 0
   yc(:,0)   = 0
   uc(:,0)   = uc(:,1)
   vc(:,0)   = vc(:,1)
   velc(:,0) = velc(:,1)

  elseif (mins .eq. 1) then

   do ct    = 0, 1
    do iang = 1, ntimes-1
     xc(ct,iang)   = gc_latlon(centre(ct,iang-1,0),centre(ct,iang-1,1),\
                               centre(ct,iang-1,0),centre(ct,iang,1),2,3)
     yc(ct,iang)   = gc_latlon(centre(ct,iang-1,0),centre(ct,iang-1,1),\
                               centre(ct,iang,0),centre(ct,iang-1,1),2,3)

   ; If storm is moving westward, zonal wind is negative
     if ( centre(ct,iang,1) .lt. centre(ct,iang-1,1) ) then
      uc(ct,iang)   = -( xc(ct,iang)/(3600) )
     else
      uc(ct,iang)   = xc(ct,iang)/(3600)
     end if

   ; If storm is moving southward, meridional wind is negative
     if ( centre(ct,iang,0) .lt. centre(ct,iang-1,0) ) then
      vc(ct,iang)   = -( yc(ct,iang)/(3600) )
     else
      vc(ct,iang)   = yc(ct,iang)/(3600)
     end if

   ; Calculate vector wind
     velc(ct,iang) = sqrt( uc(ct,iang) ^ 2.0 + vc(ct,iang) ^ 2.0)

    end do              ; End loop over times
   end do               ; End loop over storm tracker methods

 ; Set initial system velocity [0] equal to that at the second time [1]
   xc(:,0)   = 0
   yc(:,0)   = 0
   uc(:,0)   = uc(:,1)
   vc(:,0)   = vc(:,1)
   velc(:,0) = velc(:,1)

  end if

; Add metadata to new arrays
  uc@units         = "m s~S~-1~N~"
  uc@description   = "Cyclone zonal velocity"

  velc@units       = "m s~S~-1~N~"
  velc@description = "Cyclone speed"

; Create grid subsets at each time, following the storm
  if (mins .eq. 0) then

   lt_arr = new((/4,ntimes,2/),float)
   ln_arr = new((/4,ntimes,2/),float)

   lt_arr(:,:,0) = centre(:,:,0)-dist
   lt_arr(:,:,1) = centre(:,:,0)+dist
   ln_arr(:,:,0) = centre(:,:,1)-dist
   ln_arr(:,:,1) = centre(:,:,1)+dist

  elseif (mins .eq. 1) then

   lt_arr = new((/2,ntimes,2/),float)
   ln_arr = new((/2,ntimes,2/),float)

   lt_arr(:,:,0) = centre(:,:,0)-dist
   lt_arr(:,:,1) = centre(:,:,0)+dist
   ln_arr(:,:,0) = centre(:,:,1)-dist
   ln_arr(:,:,1) = centre(:,:,1)+dist

  end if

  return([/uc, vc, velc, ntimes, lt_arr, ln_arr, centre/])

end
