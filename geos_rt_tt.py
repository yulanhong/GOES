
import numpy as np

import xarray as xr

import sys

#data_dir='/data/accp/a/snesbitt/geos_tt/'

data_dir=''

datestr=sys.argv[1]

run=sys.argv[2]

fhr=sys.argv[3]

 

print('Calculating GEOS tropopause temperature: '+datestr+' Run: '+run+' Forecast hour: '+fhr)

 

m3d = xr.open_dataset('https://opendap.nccs.nasa.gov/dods/gmao/geos-cf/fcast/met_inst_1hr_g1440x721_p23/met_inst_1hr_g1440x721_p23.'+datestr+'_'+run+'z')

mds = xr.open_dataset('https://opendap.nccs.nasa.gov/dods/gmao/geos-cf/fcast/met_tavg_1hr_g1440x721_x1/met_tavg_1hr_g1440x721_x1.'+datestr+'_'+run+'z')

 

troppb = mds.isel(lev=0,time=int(fhr))['troppb']

 

t = m3d['t'].isel(time=int(fhr))

 

 

 

def interp1d_np(data, x, xi):

    return np.interp(xi, np.flip(x), np.flip(data))

 

 

#interped = interp1d_np(t.isel(lat=0,lon=0), m3d.lev,troppb.isel(lat=0,lon=0))

out = xr.apply_ufunc(

    interp1d_np,  # first the function

    t,

    m3d.lev,

    troppb/100.,  #convert from Pa to hPa to match lev units

    input_core_dims=[["lev"], ["lev"], []],

    exclude_dims=set(("lev",)),

#    output_core_dims=[["lev"]],

    vectorize=True

)

 

out.to_netcdf(data_dir+datestr+'_'+run+'_f'+fhr+'.nc')
