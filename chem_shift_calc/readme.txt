
/anaconda3/bin/python  get_chem_shift.py -h
usage: get_chem_shift.py [-h] -b BASE -l PDB_LIST -d PDB_DIR

Read in chemical shifts from nmr experiment

optional arguments:
  -h, --help            show this help message and exit
  -b BASE, --base BASE  The base of the output file.
  -l PDB_LIST, --pdb_list PDB_LIST
                        The list of pdbs to iterate through.
  -d PDB_DIR, --pdb_dir PDB_DIR
                        The directory with the pdbs to predict chemical shifts for.

python  get_chem_shift.py -l ../../../af_pred/top_pred/6VE9_mc_top.txt -d ../../../af_pred/top_pred/ -b 6VE9_mc

