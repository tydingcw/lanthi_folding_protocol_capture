#!/bin/bash

#Rosetta version
ROSETTA=/home/tydingcw/Rosetta/rosetta/source/bin/cst_info.linuxgccrelease

PDB=$1
PREFIX=$2

#$ROSETTA @cst_info.options -in:file:s ${PREFIX}*pdb* -constraints:cst_fa_file ${PDB}_noe.cst -out:file:scorefile ${PDB}_cst.sc
$ROSETTA @cst_info.options -in:file:l ${PDB}_top.txt -constraints:cst_fa_file ${PDB}_noe_neg.cst -out:file:scorefile ${PDB}_neg_cst.sc

