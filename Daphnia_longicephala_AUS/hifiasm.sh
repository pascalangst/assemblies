#!/bin/bash

#SBATCH --job-name=hifiasm                   #This is the name of your job
#SBATCH --cpus-per-task=8                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16G              #This is the memory reserved per core.
#Total memory reserved: 16GB

#SBATCH --time=120:00:00        #This is the time that your task will run
#SBATCH --qos=1week          #You will run in this queue


# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=hifiasm.out     #These are the STDOUT and STDERR files
#SBATCH --error=hifiasm.err
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
conda activate pod5

#export your required environment variables below
#################################################


#add your command lines below
#############################



hifiasm -o AUS.asm -t8 ../Pool_DE_Circular_Consensus_Sequencing_Reads/m84177_240307_192125_s3.bc2084--bc2084.hifi_reads.fastq.gz




