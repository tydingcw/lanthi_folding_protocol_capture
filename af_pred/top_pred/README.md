# Structure Prediction and Analysis for Helical Lanthipeptides

This directory contains the scropts for starting with the helical AF peptide models and introducing lanthionine rings and other non-canonical amino acids (NCAA). These structures with NCAA are predicted with Rosetta and then analyzed further and compared against the structures in the PDB.

## Scripts and Processing

- Use `lanthionine_rename.py` to rename residues in PDB files according to specific criteria

For these files, I needed to manually replace ABA with DBB: for 6VGT and 6VE9

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
Generating directories for PDB structure relax output with NMR restraints:
`echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -I % mkdir %_nmr`

Relax top predicted structures by energy with nmr constraints:
`echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -P 20 -I % bash sample_nmr_cst.sh %_free_nmr.options %`

Analyze the top peptide RMSD and calculate the NOE agreement:
```
echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -I % ~/anaconda3/bin/python get_top_nmr_cst.py -i %_free_nmr_relax.sc
echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -I % bash cst_info_nmr_opt.sh %
echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -I % bash rmsd_bb_rosetta_nmr_opt.sh %
```

MC Sampling to identify lower energy peptide conformations:
`echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -I % bash submit_mc_bluefin.sh %`

New top structures:
`echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -I % python ../../conf_gen/multi_ripp_scan/xml_test/get_quantile_list.py -i %_mc.sc -e 0.01 --outfile %_mc_top.txt`

NMR Constraint Calculation:
`echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % echo cst_info_mc.sh % %_mc_top.txt`

Plotting RMSD and energy of the MC sampled peptides - from analysis_update:
```
ls ????_mc.sc | xargs -I % echo % | sed 's/_mc.sc//' | while read -r base; do ~/anaconda3/bin/python ../gen_plots.py -i ${base}_mc.sc -c 10; done
ls ????_mc.sc | xargs -I % echo % | sed 's/_mc.sc//' | while read -r base; do bash plot_global_funnels.sh ${base}_mc_rmsd_bb.dat; done
```

Getting averaged energy and RMSD for ring:
```
echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % ~/anaconda3/bin/python /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/conf_gen/multi_ripp_scan/xml_test/calculate_avg_ring_rmsd.py -s %_mc.sc
echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % bash /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/conf_gen/multi_ripp_scan/xml_test/plot_ring_avg_funnels.sh %_mc_ring_avg_rmsd.dat
```

In conf_gen/nmr_cst/test, get chemical shift data
`~/anaconda3/bin/python  ../get_chem_shift.py -l ../../../af_pred/top_pred/6VE9_mc_top.txt -d ../../../af_pred/top_pred/ -b 6VE9_mc`

Generating Ensemble with MC refined conformations:
```
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 2KTO_noe.tbl -i 2KTO_mc_cst.sc -n 500000 -s 2KTO_mc.sc -m m -t 1.0 -e 0.01 --outfile 2KTO_mc_ensemble.txt
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 8CWX_noe.tbl -i 8CWX_mc_cst.sc -n 500000 -s 8CWX_mc.sc -m m -t 4.0 -e 0.01 --outfile 8CWX_mc_ensemble.txt
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 6VGT_noe.tbl -i 6VGT_mc_cst.sc -n 500000 -s 6VGT_mc.sc -m m -t 2.8 -e 0.01 --outfile 6VGT_mc_ensemble.txt
```

Ensemble with chemical shift
`~/anaconda3/bin/python /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/conf_gen/multi_ripp_scan/xml_test/ensemble_analysis/find_best_ensemble.py -r 6VE9_noe.tbl -i 6VE9_mc_cst.sc -n 500000 -s 6VE9_mc.sc -m m -t 8.0 --pred_shift /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/conf_gen/nmr_cst/test/6VE9_mc_shift.csv --obs_shift nmr_analysis/6VE9_shifts.txt --outfile 6VE9_mc_ensemble.txt`

RMSD BB Calculation
`echo "2KTO 6VE9 6VGT 8CWX" | sed 's/ /\n&/g' | xargs -I % bash rmsd_bb_rosetta_mc.sh %`

Pymol session generation example:
```
item=1AJ1; pymol ~/Documents/EGFR_antibodies/RiPP_design/${item}.pdb $(cat ${item}_mc_ensemble.txt | grep -v sc | tr "\n" " ");
run prep_pymol.py
```

NMR Relax Analysis
```
ls ????_free_nmr_relax.sc | xargs -I % echo % | sed 's/_free_nmr_relax.sc//' | while read -r base; do ~/anaconda3/bin/python ../gen_plots.py -i ${base}_free_nmr_relax.sc -c 10; done
ls ????_free_nmr_relax.sc | xargs -I % echo % | sed 's/_free_nmr_relax.sc//' | while read -r base; do bash plot_global_funnels.sh ${base}_free_nmr_relax_rmsd_bb.dat; done
```


