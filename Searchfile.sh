#!/bin/bash
#PBS -P Project_Name_of_Job
#PBS -q parallel8
#PBS -l select=1:ncpus=8:mpiprocs=8:mem=10GB
<<<<<<< HEAD
#PBS -e 'error.txt'
#PBS -o 'output.txt'
=======
#PBS -e '/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/error.txt'
#PBS -o '/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/output.txt'
>>>>>>> origin/master




<<<<<<< HEAD
IFS=$'\n' read -d '' -r -a lines < 'MetaData.Input.txt'
=======
IFS=$'\n' read -d '' -r -a lines < '/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/MetaData.Input.txt'
>>>>>>> origin/master

#echo "${lines[@]}"

for i in "${lines[@]}"
do
       #echo "$i"
<<<<<<< HEAD
       Rscript 'mFam-Contributions.R' "$i"
=======
       Rscript '/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/mFam-Contributions.R' "$i"
>>>>>>> origin/master
done
