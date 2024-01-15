#!/bin/bash

#SBATCH --job-name=ccsmeth                   #This is the name of your job
#SBATCH --cpus-per-task=10                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=5G              #This is the memory reserved per core.
#Total memory reserved: 16GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue
#SBATCH --partition=a100
#SBATCH --gres=gpu:2

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=ccsmeth_2.out     #These are the STDOUT and STDERR files
#SBATCH --error=ccsmeth_2.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=pascal.angst@unibas.ch        #You will be notified via email when your task ends or fails

#This job runs from the current working directory

#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


#load your required modules below
#################################
source /scicore/home/ebertd/angpas00/.bashrc 
conda activate ccsmethenv

#export your required environment variables below
#################################################


#add your command lines below
#############################

CUDA_VISIBLE_DEVICES=0 ccsmeth call_mods \
  --input m64156_230929_124402.hifi-reads.bam \
  --model_file ../model_ccsmeth_5mCpG_call_mods_attbigru2s_b21.v2.ckpt \
  --output output.hifi.call_mods \
  --threads 10 --threads_call 2 --model_type attbigru2s \
  --mode denovo




#singularity run --nv  /export/soft/singularity-containers/dorado/dorado_0.3.4.sif dorado duplex /scicore/home/ebertd/angpas00/software/dorado/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 -v pod5_pass/ > duplex_231013.bam

#/scicore/home/ebertd/angpas00/software/dorado/dorado-0.3.4-linux-x64/bin/dorado duplex /scicore/home/ebertd/angpas00/software/dorado/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 -v pod5_pass/ > duplex_230824.bam
