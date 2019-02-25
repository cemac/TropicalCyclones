;======================================================
; Function: nc_levs
; Calculate the number of vertical levels in a file
;======================================================

undef("nc_levs")
function nc_levs(fili:string)
local f, names, levs, psize

begin

; Read in all information from input file 
  f = addfile(fili,"r")

; Print names of all variables in file 
  names = getfilevarnames(f)

; Find variables that match a given string  
  zstr  = str_match_ic(names,"hybrid_ht")

; Find vertical level variable, if it exists
  if ( any (names .eq. "p") ) then 

    levs = f->p
    psize = dimsizes(levs)
    print("File contains "+psize+" pressure levels.")

  elseif ( any (names .eq. "hybrid_ht") ) then 

    levs = f->hybrid_ht
    psize = dimsizes(levs)
    print("File contains "+psize+" model levels.")

  end if 

  return([/levs,psize/])

end