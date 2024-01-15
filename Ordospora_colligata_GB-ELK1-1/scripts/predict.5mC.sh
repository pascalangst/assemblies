#!/bin/bash

#SBATCH --job-name=ccs                   #This is the name of your job
#SBATCH --cpus-per-task=64                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=2G              #This is the memory reserved per core.
#Total memory reserved: 16GB

#SBATCH --time=100:00:00        #This is the time that your task will run
#SBATCH --qos=1week          #You will run in this queue


# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=5mc.out     #These are the STDOUT and STDERR files
#SBATCH --error=5mc.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=pascal.angst@unibas.ch        #You will be notified via email when your task ends or fails

#This job runs from the current working directory

#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


#load your required modules below
#################################

#export your required environment variables below
#################################################


#add your command lines below
#############################

~/software/ccs m64156_230929_124402.subreads.bam m64156_230929_124402.hifi-reads.bam --hifi-kinetics

#~/software/jasmine m64156_230929_124402.hifi-reads.bam m64156_230929_124402.5mc.hifi-reads.bam 

