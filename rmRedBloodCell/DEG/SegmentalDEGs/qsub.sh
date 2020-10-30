#!/bin/bash
#-V:pass all environment variables to the job
#-N: job name
#-l: specify the amount of maximum memory required
#-d: working directory
work_dir=$PWD #default current working directory
echo $work_dir
#qsub -V -N qTest -l h_vmem=5G -d $work_dir test.sh
#qsub -V -N mt -l mem=10000MB -v "arg1=1000,arg2=444" -d $work_dir run1.sh
qsub -V -N qrun0  -v "arg1=0" -l mem=10000MB -d $work_dir run.sh
qsub -V -N qrun1  -v "arg1=1" -l mem=10000MB -d $work_dir run.sh
qsub -V -N qrun2  -v "arg1=2" -l mem=10000MB -d $work_dir run.sh
