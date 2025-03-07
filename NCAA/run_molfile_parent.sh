sdf=$1
charge=$2
parent=$3
tag=`basename $sdf _3D_0.sdf`

python ~/Rosetta/main/source/scripts/python/public/molfile_to_params_polymer.py \
--clobber --all-in-one-pdb --name $tag -i $sdf --partial_charges $charge --use-parent-rotamers $3
