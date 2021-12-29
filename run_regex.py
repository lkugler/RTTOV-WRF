import os
import sys
import glob
import multiprocessing as mp

import paths
RTpy = paths.python
scripts_dir = os.path.dirname(os.path.abspath(__file__))

def list_dirs(path):
    return [a for a in os.listdir(path) if os.path.isdir(os.path.join(path, a))]

def fun(filein):
    cmd = RTpy+' '+scripts_dir+'/rttov_wrf.py '+filein+' both'
    os.system(cmd) 
    print(cmd)

if __name__ == '__main__':

    regex = sys.argv[1]
    print ('got', regex)
    regex = regex.replace('\\', '')

    files = sorted(glob.glob(regex))
    if not files:
        raise ValueError('Files do not exist. Try `ls '+regex+'` to find your files!')
    pool = mp.Pool(processes=12)

    for f in files:
        print('processing:', f, flush=True)
        p = os.path.dirname(f)+'/RT_'+os.path.basename(f)+'.nc'
        if os.path.isfile(p):
            print(p, 'exists, skipping...', flush=True)
        else:
            pool.apply_async(fun, args=(f,))
    pool.close()
    pool.join()
