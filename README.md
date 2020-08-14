# RTTOV-WRF
Simulate SEVIRI satellite channels from WRF output

- Input: wrfout file (NetCDF format) written by WRF
- Output: numpy array of brightness temperatures or reflectances.


### To use RTTOV with WRF, 
1) download and compile RTTOV from [nwpsaf.eu](https://www.nwpsaf.eu/site/software/rttov/).
2) get the `pysolar` module
3) configure the python script `pyrttov_IR+VIS.py`, it sets all the assumptions for radiative transfer.


### Example call 
with path to wrfout file as argument:

`python pyrttov_IR+VIS.py path/to/wrfout_d01_2008-07-30_06:00:00`

### Output
It creates subfolders with output numpy array files in it (within the folder where the input file resides):

`./rttov_output/YYYY-mm-dd_HH:MM/RTout.YYYY-mm-dd_HH:MM.VIS06.npy`

for the channels activated in `pyrttov_IR+VIS.py`. 

### You want to contribute
Please open an issue.
Current To-Do-List:
1) change from numpy output to netcdf with xarray
2) 

### Note:
The python version calling `pyrttov_IR+VIS.py` must be the same as the one used for compiling RTTOV.
