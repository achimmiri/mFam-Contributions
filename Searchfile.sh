#!/bin/bash
#PBS -P Project_Name_of_Job
#PBS -q parallel8
#PBS -l select=1:ncpus=8:mpiprocs=8:mem=10GB
#PBS -e 'error.txt'
#PBS -o 'output.txt'




IFS=$'\n' read -d '' -r -a lines < 'MetaData.Input.txt'

#echo "${lines[@]}"

for i in "${lines[@]}"
do
       #echo "$i"
       Rscript 'mFam-Contributions.R' "$i"
done
