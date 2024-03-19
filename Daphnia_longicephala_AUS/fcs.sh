#!/bin/bash

#SBATCH --job-name=fcs_gx                   #This is the name of your job
#SBATCH --cpus-per-task=8                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 16GB

#SBATCH --time=100:00:00        #This is the time that your task will run
#SBATCH --qos=1week          #You will run in this queue
#SBATCH --partition=bigmem

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=fcs_gx.out     #These are the STDOUT and STDERR files
#SBATCH --error=fcs_gx.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=pascal.angst@unibas.ch        #You will be notified via email when your task ends or fails


#This job runs from the current working directory

#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


#load your required modules below
#################################
#source /scicore/home/ebertd/angpas00/.bashrc 
source ~/software/fcs_env/bin/activate

#export your required environment variables below
#################################################


#add your command lines below
#############################



export FCS_DEFAULT_IMAGE=~/software/ncbi-fcs/fcs-gx.sif 

python3 ~/software/ncbi-fcs/fcs.py screen genome --fasta ~/software/ncbi-fcs/fcsgx_test.fa.gz --out-dir ./gx_test_out/ --gx-db ~/software/ncbi-fcs/test-only --tax-id 6973 

python3 ~/software/ncbi-fcs/fcs.py screen genome --fasta AUS.asm.bp.p_ctg.fa --out-dir ./gx_out/ --gx-db ~/software/ncbi-fcs/gxdb --tax-id 85468 


