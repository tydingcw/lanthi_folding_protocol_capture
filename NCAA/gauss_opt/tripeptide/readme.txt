/sb/apps/corina/4.3/bin/corina_annual_x86-64_rhel6_2024_02_28.lnx -d wh,r2d DBU_Tripeptide.sdf > DBU_Tripeptide3D.sdf

Angles set according to ff19SB paper
Beta:
PyMOL>set_dihedral index 12, index 11, index 10, index 6, 135
PyMOL>set_dihedral index 11, index 10, index 6, index 3, -135
PyMOL>set_dihedral index 6, index 3, index 2, index 1, 135
PyMOL>set_dihedral index 3, index 2, index 1, index 5, -135
PyMOL>set_dihedral index 1, index 5, index 9, index 14, 135
PyMOL>set_dihedral index 5, index 9, index 14, index 15, -135

~/anaconda3/bin/python ~/scripts/gaussian_help.py -t "M052X/6-311+G(d,p)" -m gjf -i DHA_Tripeptide_beta.sdf -o DHA_Tripeptide_beta.sdf

Use gv to convert log into sdf
~/anaconda3/bin/python ~/scripts/gaussian_help.py -i DBU_Tripeptide_beta_scan.sdf -m zmat -o DBU_Tripeptide_beta_scan.gjf -d 6,3,2,1
~/miniconda3/envs/pyrosetta/bin/python ~/scripts/gaussian_rot_lib.py -l DHA_Tripeptide_beta_scan_0.log

Needed to redefine amine hydrogens based on dihedrals with atoms not in the scaned phi and psi
geom=nocrowd if necessary to bypass atoms too close

for i in {1..35}; do cp DHA_Tripeptide_beta_scan_0.gjf DHA_Tripeptide_beta_scan_${i}.gjf; sed -i "s/scan_0.chk/scan_${i}.chk/" DHA_Tripeptide_beta_scan_${i}.gjf; sed -i "s/C2 0.0/C2 ${i}0.0/" DHA_Tripeptide_beta_scan_${i}.gjf; done

#check status
ls *.log | xargs -I % sh -c 'echo % && cat % | tail -4'

