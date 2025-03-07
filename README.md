# Lanthipeptide Folding Protocol Project

This project contains various directories and scripts related to the study and analysis of lanthipeptides.

## Directory Structure

- `af_pred/`: Contains the AlphaFold predictions for four-helical lanthipeptides, as well as further analysis of these peptides. The peptides considered helical are `2KTO`, `6VE9`, `6VGT`, and `8CWX`.
- `input_pdb/`: Stores input PDB files of lanthipeptides sourced from the Protein Data Bank (PDB).
- `flex_ripp_scan/`: Includes structure prediction and analysis of non-helical lanthipeptides. The non-helical peptides are `1AJ1`, `6PQG`, `6VHJ`, `6VLJ`, `7JU9`, and `7JVF`.
- `NCAA/`: Houses the parameterization of non-canonical amino acids (NCAA) used in the study.
- `renamed_pdb/`: Contains PDB files that have been processed and renamed according to the project's conventions.
- `scripts/`: A collection of Python scripts used for various analysis and processing tasks within the project.

## Scripts

- `get_pert_rmsd_ddg.py`: Analyzes RMSD and ddG perturbations.
- `lanthionine_rename.py`: Renames residues in PDB files according to specific criteria so that they can be read in by Rosetta and lanthionine ring constraints properly applied.
- `make_nmr_cst.py`: Converts NMR restraints (tbl) to Rosetta a constraints file.
- `make_ripp_structpred_xml.py`: Creates XML files for structure prediction of RiPPs (ribosomally synthesized and post-translationally modified peptides).
- `plot_folding_funnel.py`: Plots folding funnel graphs for visualizing protein folding landscapes.

Example usage for renaming pdbs:

```bash
python ~/scripts/lanthionine_rename.py -i 7JVF.pdb -p renamed_pdb/
python ~/scripts/lanthionine_rename.py -i 6PQG.pdb -p renamed_pdb/ -f 6PQG_temp_named.pdb
