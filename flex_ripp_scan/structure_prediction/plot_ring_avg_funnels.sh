#!/bin/bash

INPUT=`readlink -e $1`
tag=`basename $INPUT .dat`

#avg_rmsd avg_energy
#-pnear
# (Å)
#
# -y_standard
/home/tydingcw/anaconda3/envs/py27/bin/python /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/af_pred/top_pred/scripts_update/plot_folding_funnel.py -input $INPUT -pnear_lambda 1.5 -pnear_kbt 0.62 -delimiter ' ' -z_str avg_energy -output ${tag}.rmsd-to-native.png -y_label "Ring Energy (REU)" -x_label "Average Ring RMSD ($\AA$)" -cmap terrain -x_str avg_rmsd -y_str avg_energy -font_size 14 -y_max 20 -exclude_axis_limits -x_standard

