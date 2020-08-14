# RTTOV-WRF
Simulate SEVIRI satellite channels from WRF output

- Input: wrfout file (NetCDF format) written by WRF
- Output: brightness temperatures or reflectances
  - numpy file
  - NetCDF file
  - xarray Dataset


### To use RTTOV with WRF, 
1) download and compile RTTOV from [nwpsaf.eu](https://www.nwpsaf.eu/site/software/rttov/).
2) get the `pysolar` module
3) configure the python script `rttov_wrf.py`, it sets all the assumptions for radiative transfer.


### Example call 
with path to wrfout file as argument:

`python rttov_wrf.py path/to/wrfout_d01 VIS xr`

### Output
It creates subfolders with output numpy array files in it (within the folder where the input file resides):

`./rttov_output/YYYY-mm-dd_HH:MM/RTout.YYYY-mm-dd_HH:MM.VIS06.npy`

for the channels activated in `rttov_wrf.py`. 

### You want to contribute
Please open an issue.

### Note:
The python version calling `rttov_wrf.py` must be the same as the one used for compiling RTTOV.
