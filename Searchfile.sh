#!/bin/bash
#PBS -P Project_Name_of_Job
#PBS -q parallel8
#PBS -l select=1:ncpus=8:mpiprocs=8:mem=10GB
#PBS -e '/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/error.txt'
#PBS -o '/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/output.txt'




IFS=$'\n' read -d '' -r -a lines < '/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/MetaData.Input.txt'

#echo "${lines[@]}"

for i in "${lines[@]}"
do
       #echo "$i"
       Rscript '/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/mFam-Contributions.R' "$i"
done
