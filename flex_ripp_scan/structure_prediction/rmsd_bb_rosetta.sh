#!/bin/bash

#Rosetta version
#ROSETTA=/home/tydingcw/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease
ROSETTA=/home/tydingcw/Rosetta/rosetta/source/bin/rosetta_scripts.default.linuxgccrelease

BASE=$1

mkdir ${BASE}_ens_rmsd

cat ${BASE}_ensemble.txt | grep -v score > ${BASE}_ensemble_process.txt

$ROSETTA -in:file:l ${BASE}_ensemble_process.txt \
-in:file:native `~/anaconda3/bin/python get_ens_min.py -e ${BASE}_ensemble.txt -s ${BASE}_scan_clean.sc`.pdb.gz \
-nstruct 1 \
-parser:protocol rmsd_bb.xml \
-use_input_sc true \
-in:file:fullatom \
-ignore_zero_occupancy false \
-linmem_ig 10 \
-out:pdb_gz \
-in:detect_disulf false \
-out:prefix ${BASE}_ens_rmsd/rmsd_bb_ \
-out:file:scorefile ${BASE}_ens_rmsd_bb.sc \
-in:file:extra_res_path /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/NCAA/uff_params/final_params/
