#!/bin/bash

INPUT=`readlink -e $1`
tag=`basename $INPUT .dat`

#avg_rmsd avg_energy
#-pnear
# (Å)
#
# -y_standard
/home/tydingcw/anaconda3/envs/py27/bin/python ~/scripts/plot_folding_funnel.py -input $INPUT -pnear_lambda 1.5 -pnear_kbt 0.62 -delimiter ' ' -z_str total_score -output ${tag}.rmsd-to-native.png -y_label "Total Energy (REU)" -x_label "Backbone RMSD ($\AA$)" -cmap terrain -x_str rmsd_bb -y_str total_score -font_size 14 -y_max 50 -exclude_axis_limits

