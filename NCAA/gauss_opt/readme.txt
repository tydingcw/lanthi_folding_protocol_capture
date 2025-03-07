Get starting conf
/sb/apps/corina/4.3/bin/corina_annual_x86-64_rhel6_2024_02_28.lnx -d wh,r2d XXX_Dipeptide.sdf > XXX_Dipeptide3D.sdf
/home/tydingcw/BCL/bcl/build/linux64_release/bin/bcl-apps-static.exe molecule:ConformerGenerator -conformation_comparer SymmetryRMSD 0.25 -max_iterations 2000 -top_models 1000 -cluster -ensemble_filenames XXX_Dipeptide3D.sdf -conformers_single_file XXX_3D_Rotamer.sdf -explicit_aromaticity
~/anaconda3/bin/python ~/scripts/min_bcl_conf.py -i XXX_3D_Rotamer.sdf

Then geom optimize
~/anaconda3/bin/python ~/scripts/gaussian_help.py -t "M052X/6-311+G(d,p)" -m gjf -i XXX_3D_Rotamer_top.sdf -o XXX_3D_Rotamer_top.sdf

