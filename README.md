# RTTOV-WRF
Simulate SEVIRI satellite channels from WRF output

- Input: wrfout file (NetCDF format) written by WRF
- Output: brightness temperatures or reflectances

### NetCDF output
with path to wrfout file as argument

Usage: `python rttov_wrf.py /path/to/wrfout (VIS|IR|both)`

E.g.: `python rttov_wrf.py /home/wrfout_d01_2008-07-30_12:00:00 both`

creates `/path/to//RT_wrfout_d01_2008-07-30_12:00:00.nc` 

Option: Use `VIS` to get VIS 6 µm reflectance, `IR` to get WV 7.3 µm & IR 10.8 µm brightness temperature or `both` to get all channels.

### xarray output 
```python

from rttov_wrf import call_pyrttov
#config = setup_IR()
config = setup_VIS()  

ds = xr.open_dataset(wrfout_path)
times = ds.Time

for t in times:
    rad = call_pyrttov(ds.sel(Time=t), config)
```

### Install
1) download and compile RTTOV from [nwpsaf.eu](https://www.nwpsaf.eu/site/software/rttov/).
2) run `python setup.py install` in the command line
3) set paths in `paths.py`
4) optional: configure the python script `rttov_wrf.py`, it sets all the assumptions for radiative transfer.


### Note:
The python version calling `rttov_wrf.py` must be the same as the one used for compiling RTTOV.
