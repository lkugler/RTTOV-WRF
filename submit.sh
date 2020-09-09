#!/bin/bash

#SBATCH -J pRTTOV

#SBATCH --account=p71386
#SBATCH --partition=mem_0384
#SBATCH --qos=p71386_0384
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --ntasks-per-core=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lukas.kugler@univie.ac.at
#SBATCH --time=30

/home/fs71386/lkugler/miniconda3/bin/python /home/fs71386/lkugler/RTTOV-WRF/loop.py

