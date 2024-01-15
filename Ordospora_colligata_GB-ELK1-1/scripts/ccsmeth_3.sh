#!/bin/bash

#SBATCH --job-name=ccsmeth                   #This is the name of your job
#SBATCH --cpus-per-task=10                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=4G              #This is the memory reserved per core.
#Total memory reserved: 16GB

#SBATCH --time=06:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue


# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=ccsmeth_3.out     #These are the STDOUT and STDERR files
#SBATCH --error=ccsmeth_3.err
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



ccsmeth align_hifi \
  --hifireads output.hifi.call_mods.modbam.bam \
  --ref GB-ELK1-1.busco.blast.self.fa \
  --output output.hifi.call_mods.modbam.pbmm2.bam \
  --threads 10


