#!/bin/bash
#-V:pass all environment variables to the job
#-N: job name
#-l: specify the amount of maximum memory required
#-d: working directory
work_dir=$PWD #default current working directory
echo $work_dir
#qsub -V -N qTest -l h_vmem=5G -d $work_dir test.sh
#qsub -V -N mt -l mem=10000MB -v "arg1=1000,arg2=444" -d $work_dir run1.sh
qsub -V -N C4 -v "arg1=4" -l mem=30000MB -d $work_dir run.sh
#qsub -V -N C2 -v "arg1=2" -l mem=30000MB -d $work_dir run.sh
#qsub -V -N C3 -v "arg1=3" -l mem=30000MB -d $work_dir run.sh
