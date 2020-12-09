#!/opt/sw/spack-0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/miniconda3-4.6.14-he6cj4ytueiygecmiasewojny57sv67s/bin/python

import os, time, warnings
import sys
import glob
from paths import path_RTTOV, path_conda

if path_conda:
    sys.path.append(path_conda)

import numpy as np
import datetime as dt
import xarray as xr
from pysolar.solar import get_altitude, get_azimuth

sys.path.append(path_RTTOV+'/wrapper')
import pyrttov

class Container(object):
    pass

def my_int(a):
    if len(a)==2 and a[0]=='0':
        return int(a[1])
    else:
        return int(a)

def expand(nprofiles, input_singlecolumn):
    return np.ones((nprofiles, 1))*input_singlecolumn

def wrftime_to_datetime(xtime):
    return dt.datetime.strptime(np.datetime_as_string(xtime, unit='s'), '%Y-%m-%dT%H:%M:%S')

def add_timezone_UTC(t):
    return dt.datetime(t.year, t.month, t.day, t.hour, t.minute, tzinfo=dt.timezone.utc)

def get_wrfout_time(ds):
    time_np = ds.XTIME.values
    time_dt = add_timezone_UTC(wrftime_to_datetime(time_np))
    return time_dt

def call_pyrttov(ds, config):
    """Run RTTOV, return xarray Dataset of reflectance or brightness temperature

    Args:
        ds (xr.Dataset): instance returned by xr.open_dataset, select one time
        config (object): see example below

    Example:
            config = setup_VIS()  # or setup_IR()
            ds_xr_out = call_pyrttov(ds_xr_in, config)
    """

    time_dt = get_wrfout_time(ds)
    print(time_dt)

    ########### input reading
    basetemp = 300.
    p = ds.PB/100.  # (ds.P + ds.PB)/100. # ignore perturbation pressure as this could break monotonicity in p, and RTTOV
    theta = ds.T+basetemp
    qv = ds.QVAPOR  #  Water vapor mixing ratio  kg kg-1
    qi = ds.QICE  # Ice mixing ratio  kg kg-1
    qc = ds.QCLOUD
    cfrac = ds.CLDFRA
    tsk = ds.TSK
    u, v = ds.U10, ds.V10
    psfc = ds.PSFC/100.
    try:  # in case that input is a wrfinput file (doesnt have these fields)
        albedo = ds.ALBEDO
        emissivity = ds.EMISS
    except AttributeError as e:
        warnings.warn(str(e)+' -> using default values!')
        albedo = 0.17 + 0*tsk
        emissivity = 0.985 + 0*tsk

    nlevels = p.shape[0]
    def reformat_profile(arr):
        # reshape to out-dims: (nprofiles, nlevels); convention from TOA to ground
        return np.flip(arr.values.reshape((nlevels, -1)).transpose(), axis=1)

    p = reformat_profile(p)
    theta = reformat_profile(theta)
    qv = reformat_profile(qv)
    qi = reformat_profile(qi)
    qc = reformat_profile(qc)
    cfrac = reformat_profile(cfrac)
    nprofiles = p.shape[0]

    # reshape to (nprofiles, nlevels)
    u = u.values.reshape((-1, 1))
    v = v.values.reshape((-1, 1))
    tsk = tsk.values.reshape((-1, 1))
    psfc = psfc.values.reshape((-1, 1))
    albedo = albedo.values.reshape((-1, 1))
    emissivity = emissivity.values.reshape((-1, 1))

    ##  conversion to Rttov
    kappa = 2/7
    tmp = theta*(p*1e-3)**kappa  # ignoring water vapor in conversion to dry temperature

    ####### actual rttov
    seviriRttov = config.seviriRttov
    irAtlas = config.irAtlas
    brdfAtlas = config.brdfAtlas

    nprofiles = p.shape[0]
    nlevels = p.shape[1]
    print('nlevels:', nlevels, 'nprofiles:', nprofiles)

    # solar angles
    lat = 45.
    lon = 0.
    sunzen = 90. - get_altitude(lat, lon, time_dt)
    sunazi = get_azimuth(lat, lon, time_dt)
    print('sunzen', sunzen, 'sunazi', sunazi)

    # surface data = lowest model half level data
    fetch = 100000.*np.ones(nprofiles)
    sfc_p_t_qv_u_v_fetch = np.stack([psfc[:,0], tmp[:,-1], qv[:,-1], u[:,0], v[:,0], fetch],
                                          axis=1).astype(np.float64)
    # sfc_u_v_fetch = np.array([[1., 1., 100000.]], dtype=np.float64)
    # sfc_u_v_fetch = expand(nprofiles, sfc_u_v_fetch)
    # sfc_p_t_qv_u_v_fetch = np.concatenate((arr[:,-1,:3], sfc_u_v_fetch), axis=-1)

    # quickfix qv
    # Q[:, :22] = 1e-6  # kg/kg  until 120 hPa
    if True:  # FIXME: is this necessary?
        qv[qv<1e-9] = 1e-9

    # print('surface p_t_qv_u_v_fetch:', sfc_p_t_qv_u_v_fetch[0,...])
    # print('lowest level  p tmp qv qc:', arr[0, -1, :])
    # print('min sp', arr[:,-1,0].min())
    # print('uppermost level p tmp qv qc:', arr[0, 0, :])

    def expand2nprofiles(n, nprof):
        # Transform 1D array to a [nprof, nlevels] array
        outp = np.empty((nprof, len(n)), dtype=n.dtype)
        for i in range(nprof):
            outp[i, :] = n[:]
        return outp

    # Declare an instance of Profiles
    myProfiles = pyrttov.Profiles(nprofiles, nlevels)

    myProfiles.GasUnits = 1  # kg/kg for mixing ratios   [ ppmv_dry=0, kg_per_kg=1, ppmv_wet=2 ]
    myProfiles.MmrCldAer = True  # kg/kg for clouds and aerosoles

    myProfiles.P = p
    myProfiles.T = tmp
    myProfiles.Q = qv
    clw = qc  # cloud liquid water
    ciw = qi  # cloud ice water

    # input on layers; (nprofiles, nlevels)
    # myProfiles.CO2 = np.ones((nprofiles,nlevels))*3.743610E+02

    # for MW channels
    # myProfiles.CLW = clw

    # quickfix
    debug = False
    if debug:
        print('cfrac', np.min(cfrac), np.max(cfrac))
    myProfiles.Cfrac = cfrac  # cloud fraction

    # WATER CLOUDS
    myProfiles.ClwScheme = 2*np.ones((nprofiles))  # "Deff"
    myProfiles.Clwde = 20*np.ones((nprofiles,nlevels))  # microns effective diameter
    # Cloud types - concentrations in kg/kg
    myProfiles.Stco = clw  # Stratus Continental STCO
    myProfiles.Stma = 0*clw  # Stratus Maritime STMA
    myProfiles.Cucc = 0*clw  # Cumulus Continental Clean CUCC
    myProfiles.Cucp = 0*clw  # Cumulus Continental Polluted CUCP
    myProfiles.Cuma = 0*clw  # Cumulus Maritime CUMA
    myProfiles.Cirr = ciw  # all ice clouds CIRR

    # icecloud[2][nprofiles]: ice_scheme, idg
    icecloud = np.array([[1, 1]], dtype=np.int32)
    myProfiles.IceCloud = expand(nprofiles, icecloud)
    myProfiles.Icede = 60 * np.ones((nprofiles,nlevels))  # microns effective diameter

    myProfiles.S2m = expand(nprofiles, sfc_p_t_qv_u_v_fetch)
    t_np = np.array([[time_dt.year,time_dt.month,time_dt.day,time_dt.hour,time_dt.minute,0]], dtype=np.int32)
    myProfiles.DateTimes = expand(nprofiles, t_np)
    # angles[4][nprofiles]: satzen, satazi, sunzen, sunazi
    angles = np.array([[45.,  180., sunzen, sunazi]], dtype=np.float64)
    myProfiles.Angles = expand(nprofiles, angles)

    # surfgeom[3][nprofiles]: lat, lon, elev
    surfgeom = np.array([[lat, lon, 0.49]], dtype=np.float64)
    myProfiles.SurfGeom = expand(nprofiles, surfgeom)
    # surftype[2][nprofiles]: surftype, watertype
    surftype = np.array([[0, 0]], dtype=np.int32)
    myProfiles.SurfType = expand(nprofiles, surftype)
    # skin[10][nprofiles]: skin T, salinity, snow_frac, foam_frac, fastem_coefsx5, specularity
    skin = np.array([[270., 35., 0., 0., 3.0, 5.0, 15.0, 0.1, 0.3, 0.]], dtype=np.float64)
    myProfiles.Skin = expand(nprofiles, skin)
    myProfiles.Skin[:,0] = tsk[:,0]

    seviriRttov.Profiles = myProfiles
    irAtlas.loadIrEmisAtlas(time_dt.month, ang_corr=True) # Include angular correction, but do not initialise for single-instrument
    brdfAtlas.loadBrdfAtlas(time_dt.month, seviriRttov) # Supply Rttov object to enable single-instrument initialisation
    brdfAtlas.IncSea = False                # Do not use BRDF atlas for sea surface types

    #######################################
    # Set up the surface emissivity/reflectance arrays and associate with the Rttov objects
    surfemisrefl_seviri = np.zeros((2, nprofiles, config.nchan), dtype=np.float64)
    surfemisrefl_seviri[0,:,:] = emissivity  # emissivity
    surfemisrefl_seviri[1,:,:] = albedo/np.pi  # reflectance
    seviriRttov.SurfEmisRefl = surfemisrefl_seviri

    albedo_emiss_by_RTTOV = False
    if albedo_emiss_by_RTTOV:
        # Surface emissivity/reflectance arrays must be initialised *before every call to RTTOV*
        # Negative values will cause RTTOV to supply emissivity/BRDF values (i.e. equivalent to
        # calcemis/calcrefl TRUE - see RTTOV user guide)
        # Call emissivity and BRDF atlases
        try:
            # Do not supply a channel list for SEVIRI: this returns emissivity/BRDF values for all
            # *loaded* channels which is what is required
            surfemisrefl_seviri[0,:,:] = irAtlas.getEmisBrdf(seviriRttov)
            surfemisrefl_seviri[1,:,:] = brdfAtlas.getEmisBrdf(seviriRttov)
        except pyrttov.RttovError as e:
            # If there was an error the emissivities/BRDFs will not have been modified so it
            # is OK to continue and call RTTOV with calcemis/calcrefl set to TRUE everywhere
            sys.stderr.write("Error calling atlas: {!s}".format(e))

    # Call the RTTOV direct model for each instrument:
    # no arguments are supplied to runDirect so all loaded channels are simulated
    try:
        t = time.time()
        seviriRttov.runDirect()
        t = time.time()-t
        print('took', int(t*10)/10., 's')
    except pyrttov.RttovError as e:
        sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))

    print('output shape:', (seviriRttov.BtRefl.shape))
    if seviriRttov.RadQuality is not None:
        print('Quality (qualityflag>0, #issues):', np.sum(seviriRttov.RadQuality > 0))


    dsout = ds.copy()
    delvars = [a for a in list(ds.variables) if not a in ['XLAT', 'XLONG']]  # 'OLR',
    dsout = dsout.drop_vars(delvars)

    # save each channel to separate file
    if debug:
        mmin = data.ravel().min()
        mmax = data.ravel().max()
        print(name, 'min:', mmin, 'max:', mmax)


    nx, ny = len(ds.west_east), len(ds.south_north)
    #data = np.zeros((ny, nx), dtype=np.float32)
    for i, name in enumerate(config.chan_seviri_names):
        data = seviriRttov.BtRefl[:,i]

        dsout[name] = (("south_north", "west_east"),
                        data.reshape(ny, nx).astype(np.float32))

    dsout.coords['time'] = time_dt.replace(tzinfo=None)  # .strftime('%Y-%m-%d_%H:%M')

    drop_lat_lon = True  # switch to False to retain coordinates from input file.
    if drop_lat_lon:
        dsout = dsout.drop_vars(['XLAT', 'XLONG'])
    return dsout


############## CONFIGURATION

def setup_IR():
    config = Container()
    seviriRttov = pyrttov.Rttov()

    chan_list_seviri = (6, 9)   #  4
    # https://nwp-saf.eumetsat.int/downloads/rtcoef_rttov12/ir_srf/rtcoef_msg_4_seviri_srf.html
    config.nchan = len(chan_list_seviri)
    config.chan_seviri_names = ('WV73', 'IR108')  # 'NIR16', 'IR39', 

    seviriRttov.FileCoef = '{}/{}'.format(path_RTTOV,
                            "/rtcoef_rttov12/rttov9pred54L/rtcoef_msg_4_seviri.dat")
    seviriRttov.FileSccld = '{}/{}'.format(path_RTTOV,
                            "/rtcoef_rttov12/cldaer_visir/sccldcoef_msg_4_seviri.dat")

    seviriRttov.Options.StoreRad = False
    seviriRttov.Options.Nthreads = 48
    seviriRttov.Options.NprofsPerCall = 840

    seviriRttov.Options.AddInterp = True
    seviriRttov.Options.AddSolar = True   # true with MFASIS
    seviriRttov.Options.AddClouds = True
    seviriRttov.Options.GridBoxAvgCloud = True
    seviriRttov.Options.UserCldOptParam = False
    seviriRttov.Options.VisScattModel = 1  # MFASIS=3  / 1 for IR sim necessary!
    seviriRttov.Options.IrScattModel  = 2
    seviriRttov.Options.OzoneData = False
    seviriRttov.Options.VerboseWrapper = False
    seviriRttov.Options.Verbose = False  # False: do not print warnings

    # ApplyRegLimits=True: Input profiles can be clipped to the regression limits when the limits are exceeded
    seviriRttov.Options.ApplyRegLimits = True

    try:
        seviriRttov.loadInst(chan_list_seviri)
    except pyrttov.RttovError as e:
        sys.stderr.write("Error loading instrument(s): {!s}".format(e))
        sys.exit(1)

    irAtlas = pyrttov.Atlas()
    irAtlas.AtlasPath = '{}/{}'.format(path_RTTOV, "/emis_data")

    brdfAtlas = pyrttov.Atlas()
    brdfAtlas.AtlasPath = '{}/{}'.format(path_RTTOV, "/brdf_data")

    config.seviriRttov = seviriRttov
    config.irAtlas = irAtlas
    config.brdfAtlas = brdfAtlas
    return config

def setup_VIS():
    config = Container()
    seviriRttov = pyrttov.Rttov()

    # select channels
    chan_list_seviri = (1,) # 2 )   # https://nwp-saf.eumetsat.int/downloads/rtcoef_rttov12/ir_srf/rtcoef_msg_4_seviri_srf.html
    config.nchan = len(chan_list_seviri)
    config.chan_seviri_names = ('VIS06', ) # 'VIS08') #, 'NIR16', 'IR39', 'WV73', 'IR108')

    # Set the options for each Rttov instance:
    # - the path to the coefficient file must always be specified
    # - turn RTTOV interpolation on (because input pressure levels differ from
    #   coefficient file levels)
    # - set the verbose_wrapper flag to true so the wrapper provides more
    #   information
    # - enable solar simulations for SEVIRI
    # - enable CO2 simulations for HIRS (the CO2 profiles are ignored for
    #   the SEVIRI and MHS simulations)

    seviriRttov.FileCoef = '{}/{}'.format(path_RTTOV,
                                          "/rtcoef_rttov12/rttov9pred54L/rtcoef_msg_4_seviri.dat")
    # CLOUD COEFFICIENTS
    seviriRttov.FileSccld = '{}/{}'.format(path_RTTOV,
                                          "/rtcoef_rttov12/cldaer_visir/sccldcoef_msg_4_seviri.dat")
    # MFASIS LOOKUPTABLE
    # seviriRttov.FileMfasisCld = '{}/{}'.format(path_RTTOV,
    #                                       "/rtcoef_rttov12/mfasis_lut/rttov_mfasis_cld_msg_4_seviri_opac.H5")
    seviriRttov.FileMfasisCld = '{}/{}'.format(path_RTTOV,
                                          "/rtcoef_rttov12/mfasis_lut/rttov_mfasis_cld_msg_4_seviri_deff.H5")

    seviriRttov.Options.StoreRad = False
    seviriRttov.Options.Nthreads = 48
    seviriRttov.Options.NprofsPerCall = 840

    seviriRttov.Options.AddInterp = True
    seviriRttov.Options.AddSolar = True
    seviriRttov.Options.AddClouds = True
    seviriRttov.Options.GridBoxAvgCloud = True
    seviriRttov.Options.UserCldOptParam = False
    seviriRttov.Options.VisScattModel = 3  # MFASIS=3
    seviriRttov.Options.IrScattModel  = 2
    seviriRttov.Options.OzoneData = False
    seviriRttov.Options.VerboseWrapper = False  #True
    seviriRttov.Options.Verbose = False  # False: do not print warnings

    # ApplyRegLimits=True: Input profiles can be clipped to the regression limits when the limits are exceeded
    seviriRttov.Options.ApplyRegLimits = True

    # Load the instruments: for HIRS and MHS do not supply a channel list and
    try:
        seviriRttov.loadInst(chan_list_seviri)
    except pyrttov.RttovError as e:
        sys.stderr.write("Error loading instrument(s): {!s}".format(e))
        sys.exit(1)

    irAtlas = pyrttov.Atlas()
    irAtlas.AtlasPath = '{}/{}'.format(path_RTTOV, "/emis_data")

    brdfAtlas = pyrttov.Atlas()
    brdfAtlas.AtlasPath = '{}/{}'.format(path_RTTOV, "/brdf_data")

    config.seviriRttov = seviriRttov
    config.irAtlas = irAtlas
    config.brdfAtlas = brdfAtlas
    return config



##########################
def calc_single_channel(channel, ds):
    """User interface to calculate sat image ad hoc

    Parameters
    ----------
    channel : str
        one out of ('VIS06', 'VIS08', 'NIR16', 'IR39', 'WV73', 'IR108')
    ds : xarray.DataArray

    Returns
    -------
    xarray.Dataset
        Contains the sat channel variable, coordinates from the input dataset.
    """
    chans = {'VIS06': 1, 'VIS08': 2, 'NIR16': 3, 'IR39': 4, 'WV73': 6, 'IR108': 9}

    ichan = chans[channel]
    if ichan >= 4:
        config = setup_IR()
    else:
        config = setup_VIS()
    config.nchan = 1
    config.chan_seviri_names = (channel, )
    dsout = call_pyrttov(ds, config)
    return dsout

if __name__ == '__main__':
    """Converts wrfout to netcdf of brightness temperature/reflectance

    Example call:
    python rttov_wrf.py /path/to/wrfout_d01 VIS

    output is one file per wrfout file, e.g. RTout_2008-07-30_18:00:00
    """
    print('>>> Usage: python rttov_wrf.py /path/to/wrfout_d01 (VIS|IR|both)')
    t0 = time.time()

    wrfout_path = sys.argv[1]
    do_both = sys.argv[2]=='both'
    do_ironly = sys.argv[2]=='IR'
    do_visonly = sys.argv[2]=='VIS'

    # do not run if output already exists
    fout = os.path.dirname(wrfout_path)+'/RT_'+os.path.basename(wrfout_path)+'.nc'
    if os.path.isfile(fout):
        print(fout, 'already exists, not running RTTOV.')
        sys.exit()
    else:
        print('running RTTOV for', wrfout_path)

    ds = xr.open_dataset(wrfout_path)
    times = ds.Time

    list_times = []
    for t in times:
        channels = []

        if do_both or do_ironly:
            config = setup_IR()
            out = call_pyrttov(ds.sel(Time=t), config)
            channels.append(out)

        if do_both or do_visonly:
            config = setup_VIS()
            out = call_pyrttov(ds.sel(Time=t), config)
            channels.append(out)

        out = xr.merge(channels)
        list_times.append(out)

    dsout = xr.concat(list_times, dim='time')
    dsout.to_netcdf(fout)
    elapsed = int(time.time() - t0)
    print(fout, 'saved, took', elapsed, 'seconds.')
