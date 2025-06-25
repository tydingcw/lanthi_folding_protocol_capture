#!/bin/bash

INPUT=`readlink -e $1`
tag=`basename $INPUT .dat`
PRE=~/scripts/
ENG=$2
VAR=$3
XLAB=$4

#echo python ${PRE}plot_folding_funnel.py -input $INPUT -pnear -pnear_lambda 1.5 -pnear_kbt 0.62 -delimiter ' ' -z_str total_score -output ${tag}.rmsd-to-best_energy.png -y_label "Energy (REU)" -cmap terrain -x_str rmsd -y_str total_score
#~/anaconda3/envs/py27/bin/python ${PRE}plot_folding_funnel.py -input $INPUT -pnear -pnear_lambda 1.5 -pnear_kbt 0.62 -delimiter ' ' -z_str total_score -output ${tag}.nmr_cst-to-best_energy.png -y_label "Energy (REU)" -cmap terrain -x_str nmr_cst -y_str total_score -x_label "NMR Constraint (REU)"
~/anaconda3/envs/py27/bin/python ${PRE}plot_folding_funnel.py -input $INPUT -delimiter ' ' -pnear -z_str $ENG -output ${tag}.rmsd-to-best_energy.png -y_label "Energy (REU)" -cmap terrain -x_str $VAR -y_str $ENG -x_label $XLAB

