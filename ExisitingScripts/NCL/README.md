# Sams NCL scripts

NCL scripts, which produce Hovmoller plots of the three-dimensional wind on model height levels [(tc_hov_uv.ncl)](tc_hov_uv.ncl) and a stamp plot of lower tropospheric windspeed on pressure levels [(tc_speed_pstamp.ncl)](tc_speed_pstamp.ncl), respectively.

# TO DO

Run and save output for comparison
* alter output
* `ncl dat=\"02T12\" opt=\"x11\" ens0=\"em11\" dist=3.0 ts=42 tf=86 mlev0=13 rmw=1 calc=1
   tclr=2 rclr=2 wclr=1 lay=1 mlev1=24 ar=1.5 nr=31 ar0=1.5 cn0=\"geo_sm\" tc_hov_uv.ncl`

## Requirements:

* `module load ncl/6.5.0/1/default`
* access to sam's user functions
* access to data
