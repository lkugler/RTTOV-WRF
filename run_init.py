import os
import sys
import glob
import datetime as dt
import multiprocessing as mp

scripts_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(scripts_dir)
import paths

def list_dirs(path):
    return [a for a in os.listdir(path) if os.path.isdir(os.path.join(path, a))]

if __name__ == '__main__':

    init = sys.argv[1]
    iens = sys.argv[2]
    force = False

    if init[-1] == '/':
        init = init[:-1]  # remove trailing /
    print ('got folder', init)

    init_time = dt.datetime.strptime(os.path.basename(init), '%Y-%m-%d_%H:%M')

    files = sorted(glob.glob(init+'/'+iens+'/wrfout_d01_20??-??-??_??:??:??'))
    format = 'wrfout_d01_%Y-%m-%d_%H:%M:%S'
    times = [dt.datetime.strptime(os.path.basename(a), format) for a in files]

    #pool = mp.Pool(processes=10)
    files_to_compute = []

    for file, time in zip(files, times):

        # save time by running only 1, 5, 10, 20, 30, 45 files
        delta_t = time - init_time
        if True: #delta_t.seconds == 90:  #  int(delta_t.seconds/60) % 10 == 0

            os.system(paths.python+' '+scripts_dir+'/rttov_wrf.py '+file+' VIS06')
    
    # pool.map_async(fun, files_to_compute)
    # pool.close()
    # pool.join()

