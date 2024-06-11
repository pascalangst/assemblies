#!/bin/bash

#SBATCH --job-name=minimap2                   #This is the name of your job
#SBATCH --cpus-per-task=12                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=3G              #This is the memory reserved per core.
#Total memory reserved: 16GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue


# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=mapping.out     #These are the STDOUT and STDERR files
#SBATCH --error=mapping.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=pascal.angst@unibas.ch        #You will be notified via email when your task ends or fails

#This job runs from the current working directory

#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


#load your required modules below
#################################
module load minimap2/2.20-GCCcore-10.3.0

#export your required environment variables below
#################################################


#add your command lines below
#############################

#minimap2 -ax map-ont -t 12 ../GiLike/RU-AST10-2.fasta duplex_230824.fq > KZ313cfsimGilike.similis.ont.mapped.sam



#minimap2 -ax map-ont -t 12 IL-YERU-16.mask.fa duplex_230824.fq > KZ313cfsimGilike.ILYERU.ont.mapped.sam


minimap2 -ax map-hifi -t 12 ../BEOM33/29082016_Xinb3_ref.fasta ../Pool_RU-NE1-34_KZ-27-1_KZ-23-1/DE3_pool_Sept_Circular_Consensus_Sequencing_Reads/m64156_230928_014710.bc2006--bc2006.hifi_reads.fastq.gz > KZ-23-1_Xinb3.sam



module purge
module load SAMtools/1.7-goolf-1.7.20

samtools flagstat KZ-23-1_Xinb3.sam > KZ-23-1_Xinb3.sam.flagstat

samtools view -Sb -f 4 KZ-23-1_Xinb3.sam | \
samtools sort -n -o KZ-23-1_Xinb3.hifi.unmapped.bam -
samtools fastq KZ-23-1_Xinb3.hifi.unmapped.bam -0 KZ-23-1_Xinb3.hifi.unmapped.fastq.gz -n

rm KZ-23-1_Xinb3.sam





