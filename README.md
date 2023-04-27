# RTTOV-WRF
Simulate SEVIRI satellite channels from WRF output

- Input: wrfout file (NetCDF format) written by WRF
- Output: brightness temperatures or reflectances

## Usage examples
### Use on command line: NetCDF output

Usage: `python rttov_wrf.py /path/to/wrfout (VIS|IR|both)`

For example, `python rttov_wrf.py /somepath/wrfout_d01_2008-07-30_12:00:00 both` 

creates `/somepath/RT_wrfout_d01_2008-07-30_12:00:00.nc` 

Option: Use `VIS` to get VIS 0.6 µm reflectance, `IR` to get WV 6.2 µm, 7.3 µm and IR 10.8 µm brightness temperature or `both` to get all channels. 
You can add or remove channels in `rttov_wrf.py`.

### Use within python: Xarray output 
```python

from rttov_wrf import call_pyrttov
#config = setup_IR()
config = setup_VIS()  

ds = xr.open_dataset(wrfout_path)
times = ds.Time

for t in times:
    rad = call_pyrttov(ds.sel(Time=t), config)
```
### Run all ensemble members wrfout files from one forecast init (12z) at 13z with:

`python run_pattern.py /path_to_exp/2008-07-30_12\:00/*/wrfout_d01_2008-07-30_13\:00\:00`


## Install
1) Download and compile RTTOV from [nwpsaf.eu](https://www.nwpsaf.eu/site/software/rttov/).
2) Download RTTOV-WRF for example by running `git clone https://github.com/lkugler/RTTOV-WRF.git`
3) Run `cd RTTOV-WRF; pip install -e .` in the command line
4) Set paths in `paths.py`
5) optional: configure the python script `rttov_wrf.py`, it sets various assumptions for radiative transfer.
6) when running `rttov_wrf.py`, ensure that you have loaded the same libraries which you used to install RTTOV. For me, this is `intel-mpi intel netcdf netcdf-fortran zlib hdf5` (on VSC: `module purge; module load intel-mpi/2019.3 intel/19.1.0 netcdf/4.7.0-intel-19.0.5.281-75t52g6 netcdf-fortran/4.4.5-intel-19.0.5.281-qye4cqn zlib/1.2.11-intel-19.1.0.166-hs6m2qh  hdf5/1.10.5-intel-19.0.5.281-qyzojtm`)

In order to run, it needs:
1) compiled RTTOV / pyrttov (the python wrapper); with the path of pyrttov in PYTHONPATH (is ensured in `rttov_wrf.py`)
2) installed RTTOV-WRF (including dependencies)
3) loaded libraries for RTTOV

### Dependencies
Will be installed by the pip command.

### Note:
If you receive an `ImportError` similar to this one
```
ImportError: libpython3.7m.so.1.0: cannot open shared object file: No such file or directory
```
at line `import pyrttov`, then your python can not find the pyrttov module. You must call `rttov_wrf.py` with the same python as was used for compiling RTTOV.

### Documentation
For the configuration of RTTOV or in case there are errors with RTTOV, consult the documentation of RTTOV or the RTTOV-python-wrapper.
(https://nwp-saf.eumetsat.int/site/software/rttov/documentation/)
