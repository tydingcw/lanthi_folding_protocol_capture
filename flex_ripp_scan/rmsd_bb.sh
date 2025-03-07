#!/bin/bash

#Rosetta version
#ROSETTA=/home/tydingcw/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease
ROSETTA=/home/tydingcw/Rosetta/rosetta/source/bin/rosetta_scripts.default.linuxgccrelease

BASE=$1

mkdir ${BASE}_native

$ROSETTA -in:file:l ${BASE}_native_list.txt \
-in:file:native /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/renamed_pdb/${BASE}_1_rename.pdb \
-nstruct 1 \
-parser:protocol rmsd_bb.xml \
-use_input_sc true \
-in:file:fullatom \
-ignore_zero_occupancy false \
-linmem_ig 10 \
-out:pdb_gz \
-in:detect_disulf false \
-out:prefix ${BASE}_native/rmsd_bb_ \
-in:file:extra_res_path /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/NCAA/uff_params/final_params/
