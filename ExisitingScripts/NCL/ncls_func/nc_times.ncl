;====================================================
; Function: nc_times
; Calculate the number of analysis times in a file
;====================================================

undef("nc_times")
function nc_times(fili:string)
local f, names, time, tsize

begin

  ntimes = 120			       ; Number of times 

; Read in all information from input file 
  f = addfile(fili,"r")

; Print names of all variables in file 
  names = getfilevarnames(f)

; Find 'time' variable
  if ( any (names .eq. "t") ) then 
    time = f->t
    tsize = dimsizes(time)
;    print("File contains "+tsize+" analysis times.")
  end if 

  return([/tsize,time/])

end