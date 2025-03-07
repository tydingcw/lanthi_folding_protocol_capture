#!/bin/bash

#Rosetta version
ROSETTA=/home/tydingcw/Rosetta/rosetta/source/bin/cst_info.linuxgccrelease

PDB=$1
PREFIX=/home/tydingcw/Documents/EGFR_antibodies/RiPP_design/renamed_pdb/

$ROSETTA @cst_info.options -in:file:s ${PREFIX}*${PDB}*pdb -constraints:cst_fa_file ${PDB}_noe.cst -out:file:scorefile ${PDB}_cst_nmr.sc

