for item in *.sdf; do /sb/apps/corina/4.3/bin/corina_annual_x86-64_rhel6_2024_02_28.lnx -d wh,r2d $item > `basename $item .sdf`_3D.sdf; done
