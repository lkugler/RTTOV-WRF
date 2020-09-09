#!/opt/sw/spack-0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/miniconda3-4.6.14-he6cj4ytueiygecmiasewojny57sv67s/bin/python

import os, time
import sys
import glob
import multiprocessing as mp

RTpy = '/opt/sw/spack-0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/miniconda3-4.6.14-he6cj4ytueiygecmiasewojny57sv67s/bin/python'

def fun(filein):
    fout = filein.replace('wrfout_d01_', 'RTout.')+'.nc'
    if os.path.isfile(fout):
        print(fout, 'exists, skip RTTOV.')
    else:
        print('running', filein)
        os.system(RTpy+' ~/RTTOV-WRF/rttov_wrf.py '+filein+' both xr')

if __name__ == '__main__':
    datadir = '/gpfs/data/fs71386/lkugler/sim_archive/'
    exps = ['exp_v1.11_LMU_filter2/',
            ]

    nature = 'exp_v1.11_LMU_nature/2008-07-30_06:00/2'
    n_ens = 40

    for exp in exps:
        inits = sorted(glob.glob(datadir+exp+'/20*'))
        print('inits', inits)
        pool = mp.Pool(processes=8)  # 8 * 6 cores = 48

        # nature run
        naturedir = os.path.join(datadir, nature)
        wrfout_files = sorted(glob.glob(naturedir+'/wrfout_d01_20*'))
        for f in wrfout_files:
            pool.apply_async(fun, args=(f,))
        pool.close()
        pool.join()


        for initdir in inits:
            for iens in range(1, n_ens+1):
                wrfout_files = sorted(glob.glob(initdir+'/'+str(iens)+'/wrfout_d01_20*'))
                print('wrfout_files')
                for f in wrfout_files:
                    pool.apply_async(fun, args=(f,))
        pool.close()
        pool.join()
