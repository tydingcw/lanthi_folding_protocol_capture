#!/bin/sh
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1G
#SBATCH --exclude bf[49-51]
#SBATCH --time=0
#SBATCH --output=slurm/%A_%a.out
#SBATCH --job-name="ripp_folding"
#SBATCH --array=1-11000

#${SLURM_ARRAY_TASK_ID}
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

#export LD_LIBRARY_PATH="/dors/meilerlab/apps/Linux2/x86_64/lib64:/dors/meilerlab/apps/Linux2/x86_64/lib:/dors/meilerlab/home/tydingcw/Rosetta/main/source/build/external/release/linux/3.10/64/x86/gcc/5.2/bcl/"
#export LD_LIBRARY_PATH=/sb/meilerapps/Linux2/x86_64/gcc/9.1.0/lib64:/sb/meilerapps/Linux2/x86_64/gcc/9.1.0/lib:$LD_LIBRARY_PATH

#Rosetta version
#ROSETTA=/home/tydingcw/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease
#ROSETTA=/home/tydingcw/Rosetta/main/source/build/src/release/linux/3.10/64/x86/gcc/5.2/default/rosetta_scripts.default.linuxgccrelease
#ROSETTA=/home/tydingcw/Rosetta/rosetta/source/bin/rosetta_scripts.default.linuxgccrelease
ROSETTA=/home/tydingcw/Rosetta/rosetta_blue/source/bin/rosetta_scripts.default.linuxgccrelease

#OPTIONS=$1
PRE=$1
SUFFIX=$2
OPTIONS=${PRE}${SLURM_ARRAY_TASK_ID}${SUFFIX}
#xml_file=test_6vlj.xml
#CSV=quick_scan.txt
#CSV=sample_instruct2.csv
#XML=predict_nmr_gen_6vlj.xml
#PDB=../renamed_pdb/6VLJ_1_rename.pdb
#PREFIX=sample_instruct_struct/
#PREFIX=sample_struct/

#tag=`basename $PDB .pdb`_rename.pdb
#echo $XML $REC_LEN $RIPP_LEN

$ROSETTA @${OPTIONS}

