#!/bin/bash

#Rosetta version
#ROSETTA=/home/tydingcw/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease
#ROSETTA=/home/tydingcw/Rosetta/main/source/build/src/release/linux/3.10/64/x86/gcc/5.2/default/rosetta_scripts.default.linuxgccrelease
#ROSETTA=/home/tydingcw/Rosetta/main/source/build/src/debug/linux/3.10/64/x86/gcc/5.2/default//rosetta_scripts.default.linuxgccdebug
ROSETTA=/home/tydingcw/Rosetta/rosetta/source/bin/rosetta_scripts.default.linuxgccrelease
#ROSETTA=/home/tydingcw/Rosetta/rosetta/source/bin/rosetta_scripts.default.linuxgccdebug

OPTIONS=$1
BASE=$2
#xml_file=test_6vlj.xml
#CSV=quick_scan.txt
#CSV=sample_instruct2.csv
#XML=predict_nmr_gen_6vlj.xml
#PDB=../renamed_pdb/6VLJ_1_rename.pdb
#PREFIX=sample_instruct_struct/
#PREFIX=sample_struct/

#tag=`basename $PDB .pdb`_rename.pdb
#echo $XML $REC_LEN $RIPP_LEN

dir=`basename $OPTIONS .options`
mkdir -p $dir

#$ROSETTA -in:file:s /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/renamed_pdb/${BASE}_*.pdb @${OPTIONS}
$ROSETTA -in:file:s ${BASE} @${OPTIONS}
#$ROSETTA -in:file:l ${BASE}_top.txt @${OPTIONS}

