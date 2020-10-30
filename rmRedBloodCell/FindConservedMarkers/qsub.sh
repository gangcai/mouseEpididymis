#!/bin/bash
#-V:pass all environment variables to the job
#-N: job name
#-l: specify the amount of maximum memory required
#-d: working directory
work_dir=$PWD #default current working directory
echo $work_dir
#qsub -V -N qTest -l h_vmem=5G -d $work_dir test.sh
#qsub -V -N mt -l mem=10000MB -v "arg1=1000,arg2=444" -d $work_dir run1.sh
for arg in {0..7}
do

	echo $arg
	qsub -V -N fcons$arg -v "arg1=$arg" -l mem=10000MB -d $work_dir run.sh
done
