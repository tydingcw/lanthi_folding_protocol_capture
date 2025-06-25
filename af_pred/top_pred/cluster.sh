#!/bin/bash

#ROSETTA=/dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin/energy_based_clustering.default.linuxgccrelease
ROSETTA=/home/tydingcw/Rosetta/main/source/bin/energy_based_clustering.default.linuxgccrelease
#/home/tydingcw/Rosetta/main/source/build/src/release/linux/3.10/64/x86/gcc/5.2/default/energy_based_clustering.default.linuxgccrelease

PDB=$1
ENS=$2
PREFIX=$3
OUTDIR=$4

mkdir $OUTDIR

#$ROSETTA @ /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/af_pred/top_pred/cluster.options -in:file:s ${PREFIX}*pdb* -out:prefix $OUTDIR/${PDB} -cluster:energy_based_clustering:alternative_score_file ${PDB}_scores.txt
$ROSETTA @ /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/af_pred/top_pred/cluster.options -in:file:l $ENS -out:prefix $OUTDIR/${PDB} -cluster:energy_based_clustering:alternative_score_file ${PDB}_scores.txt
