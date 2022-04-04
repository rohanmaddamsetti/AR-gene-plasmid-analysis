#!/bin/bash
#SBATCH --job-name=run_ICRA       # create a short name for your job
#SBATCH --output=slurm-%A.%a.out  # stdout file
#SBATCH --error=slurm-%A.%a.err   # stderr file
#SBATCH --nodes=1                 # node count
#SBATCH --ntasks=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=1         # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=32G         # memory per cpu-core (need at least 20G)
#SBATCH --time=24:00:00           # total run time limit (HH:MM:SS)
#SBATCH --array=0-239             # job array with index values 0-239 (240 total)
#SBATCH --partition=scavenger     # run in the scavenger queue.
#SBATCH --mail-type=FAIL,END      # can send email on job start, end, fail, or all.
#SBATCH --mail-user=rohan.maddamsetti@duke.edu

## run_ICRA_on_DIABIMMUNE_jobarray.sh by Rohan Maddamsetti.
## I used the template jobarray script from here:
## https://researchcomputing.princeton.edu/support/knowledge-base/slurm#arrays
## Usage: sbatch run_ICRA_on_DIABIMMUNE_jobarray.slurm

## IMPORTANT: before running this script, conda activate sgvfinder-env.


echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

python run_ICRA_on_DIABIMMUNE.py
