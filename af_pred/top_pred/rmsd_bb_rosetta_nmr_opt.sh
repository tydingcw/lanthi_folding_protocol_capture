#!/bin/bash

#Rosetta version
#ROSETTA=/home/tydingcw/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease
ROSETTA=/home/tydingcw/Rosetta/rosetta/source/bin/rosetta_scripts.default.linuxgccrelease

BASE=$1

mkdir -p ${BASE}_ens_rmsd

#cat ${BASE}_ensemble.txt | grep -v score > ${BASE}_ensemble_process.txt

$ROSETTA -in:file:l ${BASE}_nmr_opt.txt \
-in:file:native `cat ${BASE}_nmr_opt.txt | head -n 1` \
-nstruct 1 \
-parser:protocol rmsd_bb.xml \
-use_input_sc true \
-in:file:fullatom \
-ignore_zero_occupancy false \
-linmem_ig 10 \
-in:detect_disulf false \
-out:file:score_only ${BASE}_ens_nmr_opt_rmsd_bb.sc \
-in:file:extra_res_path /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/NCAA/uff_params/final_params/

#-out:pdb_gz \
#-out:prefix ${BASE}_ens_rmsd/rmsd_bb_ \

