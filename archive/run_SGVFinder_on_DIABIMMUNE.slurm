#!/bin/bash
#SBATCH --job-name=run_SGVFinder  # create a short name for your job
#SBATCH --output=slurm-%A.%a.out  # stdout file
#SBATCH --error=slurm-%A.%a.err   # stderr file
#SBATCH --nodes=1                 # node count
#SBATCH --ntasks=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=240       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=8G          # memory per cpu-core
#SBATCH --time=24:00:00           # total run time limit (HH:MM:SS)
#SBATCH --partition=scavenger     # run in the scavenger queue.
#SBATCH --mail-type=FAIL,END      # can send email on job start, end, fail, or all.
#SBATCH --mail-user=rohan.maddamsetti@duke.edu

## run_SGVFinder_on_DIABIMMUNE.slurm by Rohan Maddamsetti.
## I used the template jobarray script from here:
## https://researchcomputing.princeton.edu/support/knowledge-base/slurm#arrays
## Usage: sbatch run_SGVFinder_on_DIABIMMUNE.slurm

## IMPORTANT: before running this script, conda activate sgvfinder-env.

echo "Executing on the machine:" $(hostname)

python run_SGVFinder_on_DIABIMMUNE.py
