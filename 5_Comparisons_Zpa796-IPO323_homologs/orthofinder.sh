#!/bin/bash

#BATCH --job-name=orthofinder #Give your job a name.

#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=25 #Multithreading.
#SBATCH --time=96:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=100G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --partition=highmem #Request a specific partition for the resource allocation.


##########################
# Module load:
#module load $MODULE_NAME 

#########################
# Run commands from here:


orthofinder -f <path2/folder_with_proteomes> -t 25 -M msa -S diamond






