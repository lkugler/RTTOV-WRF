#!/bin/bash
#PBS -N logpbs
#PBS -A NMMM0063
#PBS -j oe
#PBS -q develop
#PBS -l walltime=00:05:00
#PBS -l select=1:ncpus=128:mpiprocs=1:mem=230GB

module load conda 
eval "$(mamba shell hook --shell bash)"
mamba activate lkugler-Py310

cd /glade/u/home/lkugler/RTTOV-WRF/
python rttov_mpas.py /glade/work/swei/projects/hydrosat/tests/test_pmo/mpasout.nc /glade/work/lkugler/hydrosat/4p16s_3km/obs.2024050111.nc /glade/work/lkugler/hydrosat/RTout/4p16s_3km_obs.2024050111.nc