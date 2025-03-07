#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=p_csb_meiler
#SBATCH --cpus-per-task=16
#SBATCH --mem=9G
#SBATCH --time=4-00:00:00
#SBATCH --array=0-35
#SBATCH --job-name=gaussian_pes
#SBATCH --output=slurm/pes_%A_%a.out

#FILE=$1

job=`echo DBU_Tripeptide_beta_scan_${SLURM_ARRAY_TASK_ID}.gjf`

#echo g16 disulfide_SS_scan_opt_step${SLURM_ARRAY_TASK_ID}.gjf
echo g16 $job
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load Gaussian/16.B.01
#g16 disulfide_SS_scan_opt_step${SLURM_ARRAY_TASK_ID}.gjf
g16 $job
