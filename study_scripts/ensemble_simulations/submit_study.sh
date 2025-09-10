#!/bin/bash
#SBATCH --qos=short
#SBATCH --job-name=run-YOUR_JOB
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --chdir=YOUR_PATH

echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "$SLURM_CPUS_PER_TASK cpus/task"
echo "------------------------------------------------------------"

module load julia/1.10.0

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun julia YOUR_SCRIPT.jl 
