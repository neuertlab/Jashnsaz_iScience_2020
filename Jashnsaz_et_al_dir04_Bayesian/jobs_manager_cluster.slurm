#!/bin/bash
#SBATCH --mail-user=jashnsh@vanderbilt.edu
#SBATCH --mail-type=ALL

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --mem=16G
#SBATCH --array=1-30
#SBATCH --output=Out_%A_%a.txt

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load MATLAB
matlab -nodisplay -nosplash -r "run_Bayes(${SLURM_ARRAY_TASK_ID});quit;"