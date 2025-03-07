# Structure Prediction and Analysis for Helical Lanthipeptides

This directory contains the scropts for starting with the helical AF peptide models and introducing lanthionine rings and other non-canonical amino acids (NCAA). These structures with NCAA are predicted with Rosetta and then analyzed further and compared against the structures in the PDB.

## Scripts and Processing

- Use `lanthionine_rename.py` to rename residues in PDB files according to specific criteria

For these files, I needed to manually replace ABA with DBB for 6VGT and 6VE9

To generate the NMR constraint files for Rosetta:
`for name in 2KTO 8CWX 6VE9 6VGT; do ~/anaconda3/bin/python ~/scripts/make_nmr_cst.py -i ${name}_noe.tbl -s ${name}_1_rename_pep.txt -c ../../renamed_pdb/${name}_conn.txt; done`

To pick fragments with the Rosetta fragment picker:
```
Making Fragments:
/sb/meilerapps/scripts/rosetta_tools/make_fragments.pl pepwt.fasta 
for file in ?????.jufo_ss; do echo -e "# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)\n\n$(cat "$file")" > temp && mv temp "$file"; done
ls ?????.fasta | xargs -I % basename % .fasta | xargs -I % bash fragment_picker.command.sh %
```
For fragment sampling, we adopt the score function from the ConfChangeMover

Commands to run the relaxation of PDB structures and structure prediction sampling:
```
ls ????_f???.options ????_bb.options | xargs -n 1 -P 12 -I % bash -c 'bash sample_native.sh % `echo % | cut -c-4`'
ls ????_af.options | xargs -n 1 -P 4 -I % bash sample_instruct.sh %
ls ????_samp.options | xargs -n 1 -P 4 -I % bash sample_instruct.sh %
```

Getting the NMR restraint agreement for each NOE term for the PDB structures:
`ls ????.span | xargs -I % echo % | sed 's/.span//' | while read -r base; do bash cst_info_nmr.sh $base; done`

Getting the NMR restraint agreement for each NOE term for every predicted structure:
`echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % echo cst_info.sh % %/`

Generating plots (from the analysis directory):
```
ls ????_af_sample.sc | xargs -I % echo % | sed 's/_af_sample.sc//' | while read -r base; do ~/anaconda3/bin/python ../gen_plots.py -i ${base}_af_sample.sc -c 10; done
ls ????_af_sample.sc | xargs -I % echo % | sed 's/_af_sample.sc//' | while read -r base; do bash plot_global_funnels.sh ${base}_af_sample_rmsd_bb.dat; done
echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % python calculate_avg_ring_rmsd.py -s %_af_sample.sc
echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % bash plot_ring_avg_funnels.sh %_af_sample_ring_avg_rmsd.dat
```

Generating Ensembles:
```
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 2KTO_noe.tbl -i 2KTO_cst.sc -n 500000 -s 2KTO_af_sample.sc -m m -t 0.8
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 8CWX_noe.tbl -i 8CWX_cst.sc -n 500000 -s 8CWX_af_sample.sc -m m -t 4.0
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 6VGT_noe.tbl -i 6VGT_cst.sc -n 500000 -s 6VGT_af_sample.sc -m m -t 3.0
```

Chemical Shift Prediciton in nmr_analysis/chem_shift - I used SPARTA+ for the chemical shift calculations

Generating Ensemble with Chemical Shift:
`python ensemble_analysis/find_best_ensemble.py -r 6VE9_noe.tbl -i 6VE9_cst.sc -n 500000 -s 6VE9_af_sample.sc -m m -t 5.0 --pred_shift nmr_analysis/chem_shift/6VE9_shift.csv --obs_shift nmr_analysis/6VE9_shifts.txt`

RMSD calculations:
```
ls ????.cst | xargs -I % echo % | sed 's/.cst//' | while read -r base; do ls ../../renamed_pdb/${base}_*pdb > "${base}_native_list.txt"; done
ls ????.span | xargs -I % echo % | sed 's/.span//' | while read -r base; do bash rmsd_bb.sh $base; done
ls ????.span | xargs -I % echo % | sed 's/.span//' | while read -r base; do bash rmsd_bb_rosetta.sh $base; done
```

