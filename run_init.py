import os
import sys
import glob
import datetime as dt
import multiprocessing as mp

scripts_dir = os.path.dirname(os.path.abspath(__file__))

def list_dirs(path):
    return [a for a in os.listdir(path) if os.path.isdir(os.path.join(path, a))]

def fun(filein):
    cmd = 'python '+scripts_dir+'/rttov_wrf.py '+filein+' both'
    os.system(cmd) 
    print(cmd)


if __name__ == '__main__':

    init = sys.argv[1]

    if init[-1] == '/':
        init = init[:-1]  # remove trailing /
    print ('got folder', init)

    init_time = dt.datetime.strptime(os.path.basename(init), '%Y-%m-%d_%H:%M')

    files = sorted(glob.glob(init+'/*/wrfout_d01_20??-??-??_??:??:00'))
    format = 'wrfout_d01_%Y-%m-%d_%H:%M:00'
    times = [dt.datetime.strptime(os.path.basename(a), format) for a in files]

    pool = mp.Pool(processes=10)

    for file, time in zip(files, times):

        # save time by running only 1, 5, 10, 20, 30, 45 files
        delta_t = time - init_time
        
        #if int(delta_t.seconds/60) in [1, 5, 10, 20, 30, 45, 60, 75]:

        #if int(delta_t.seconds/60) % 10 == 0 or int(delta_t.seconds/60) in [1, 5]:

        #if ((int(time.minute) in [5, 10, 15, 20, 30, 35, 40, 45, 50])
        #   or (int(time.minute) == 0 and int(delta_t.seconds) > 0)):

        if True:

            # save time, dont rerun existing files
            p = os.path.dirname(file)+'/RT_'+os.path.basename(file)+'.nc'
            if not os.path.isfile(p):
                print(p, 'does not yet exist', flush=True)
                pool.apply_async(fun, args=(file,))
    pool.close()
    pool.join()

