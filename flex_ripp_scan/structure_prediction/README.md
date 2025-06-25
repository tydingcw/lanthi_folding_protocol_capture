# Structure Prediction Analysis

This directory contains scripts and data for the analysis of structure prediction results for lanthipeptides.

- Generate the initial XML and options for structure prediction using `make_ripp_structpred_xml.py`:

  ```bash
  ~/anaconda3/bin/python ..//make_ripp_structpred_xml.py
  ```

## Analysis

Generate plots for each `.sc` score file:

```bash
for item in *.sc; do
    ~/anaconda3/bin/python ~/scripts/gen_plots.py -i $item
done
```

Finalize plots from the analysis directory:

```bash
bash ../pnear.sh 6VHJ_scan_ring1_rmsd_bb.dat ring1_eng ring1_rmsd_bb Ring1/ RMSD/ Backbone
```

Get NMR restraint values:

```bash
echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | tr ' ' '\n' | while read -r dir; do
    find "$dir/" -name "*pdb*" > "$dir.txt"
done
echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % echo cst_info.sh % %/
```

Calculate negative constraints:

```bash
~/anaconda3/bin/python calculate_negative_cst.py -p 6VLJ/1_init_6VLJ_1_rename_0005.pdb.gz -c 6VLJ_noe.cst
# Repeat for other peptides as needed
```

Clean up score files:

```bash
for file in ????_scan.sc; do
    ~/anaconda3/bin/python clean_scorefile.py -s $file
done
```

Get the structures in the top quantile:

```bash
~/anaconda3/bin/python get_quantile_list.py -i 1AJ1_scan_clean.sc -e 0.025
# Output: Energy threshold is: 59.695575000000005
# Repeat for other peptides as needed
```

Get the negative constraints for the top quantile:

```bash
echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % echo cst_info_neg.sh % %/
```

Ensemble analysis:

```bash
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 1AJ1_noe.tbl -i 1AJ1_cst.sc -n 500000 -s 1AJ1_scan_clean.sc -m m -e 0.025 -t 8.0
# Repeat for other peptides as needed
```

NMR constraint information:

```bash
ls ????.cst | xargs -I % echo % | sed 's/.cst//' | while read -r base; do
    bash cst_info_nmr.sh $base
done
```
Get averaged energy and RMSD for the ring:

```bash
echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % ~/anaconda3/bin/python calculate_avg_ring_rmsd.py -s %_scan_clean.sc
echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % bash plot_ring_avg_funnels.sh %_scan_clean_ring_avg_rmsd.dat
```

Global RMSD analysis:

```bash
echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % bash ../plot_global_funnels.sh %_scan_clean_rmsd_bb.dat
echo "1AJ1 6PQG" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % bash ../plot_global_funnels_yraise.sh %_scan_clean_rmsd_bb.dat
```

Filter by atropisomers
```
~/anaconda3/bin/python ensemble_analysis/atripisomer_parse.py -p 1AJ1_top.txt -d 17_CA-19_CA-15_CA-9_CA 9_CA-11_CA-7_CA-17_CA -r -f neg_ens/1AJ1_ensemble.txt -o 1AJ1_ensemble_filt.csv
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 1AJ1_noe.tbl -i 1AJ1_cst.sc -n 500000 -s 1AJ1_scan_clean.sc -m m -e 0.025 -t 8.0 -f 1AJ1_ensemble_filt.csv
~/anaconda3/bin/python ensemble_analysis/atripisomer_parse.py -p 6PQG_top.txt -d 16_CA-17_CA-14_CA-11_CA 15_CA-17_CA-13_CA-10_CA -r -f neg_ens/6PQG_ensemble.txt -o 6PQG_ensemble_filt.csv
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 6PQG_noe.tbl -i 6PQG_cst.sc -n 500000 -s 6PQG_scan_clean.sc -m m -e 0.025 -t 8.0 -f 6PQG_ensemble_filt.csv
```

Generating directories for PDB structure relax output with NMR restraints:
`echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -I % mkdir %_nmr`

Relax top predicted structures by energy with nmr constraints:
`echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -P 20 -I % bash sample_nmr_cst.sh %_free_nmr.options %`

Analyze the top peptide RMSD and calculate the NOE agreement:
```
echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -I % ~/anaconda3/bin/python get_top_nmr_cst.py -i %_free_nmr_relax.sc
echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -I % bash cst_info_nmr_opt.sh %
echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -I % bash rmsd_bb_rosetta_nmr_opt.sh %
```

MC Sampling to identify lower energy peptide conformations:
`echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -I % bash submit_mc_bluefin.sh %`

New top structures:
`echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -I % ~/anaconda3/bin/python get_quantile_list.py -i %_mc.sc -e 0.025 --outfile %_mc_top.txt`

NMR Constraint Calculation:
`echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % echo cst_info_mc.sh % %_mc_top.txt`

Plotting RMSD and energy of the MC sampled peptides - from analysis_update:
```
ls ????_free_nmr_relax.sc | xargs -I % echo % | sed 's/_free_nmr_relax.sc//' | while read -r base; do ~/anaconda3/bin/python ../gen_plots.py -i ${base}_free_nmr_relax.sc -c 10; done
ls 1AJ1_free_nmr_relax.sc 6PQG_free_nmr_relax.sc | xargs -I % echo % | sed 's/_free_nmr_relax.sc//' | while read -r base; do ~/anaconda3/bin/python ../gen_plots.py -i ${base}_free_nmr_relax.sc -c 30; done
ls ????_free_nmr_relax.sc | xargs -I % echo % | sed 's/_free_nmr_relax.sc//' | while read -r base; do bash plot_global_funnels.sh ${base}_free_nmr_relax_rmsd_bb.dat; done
```

Getting averaged energy and rmsd for ring:
```
echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % ~/anaconda3/bin/python /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/conf_gen/multi_ripp_scan/xml_test/calculate_avg_ring_rmsd.py -s %_mc.sc
echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -n1 -P 6 -I % bash /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/conf_gen/multi_ripp_scan/xml_test/plot_ring_avg_funnels.sh %_mc_ring_avg_rmsd.dat
```

Generating Ensemble with MC refined conformations:
```
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 6VHJ_noe.tbl -i 6VHJ_mc_cst.sc -n 500000 -s 6VHJ_mc.sc -m m -e 0.025 -t 22.0 --outfile 6VHJ_mc_ensemble.txt
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 6VLJ_noe.tbl -i 6VLJ_mc_cst.sc -n 500000 -s 6VLJ_mc.sc -m m -e 0.025 -t 21.0 --outfile 6VLJ_mc_ensemble.txt
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 7JU9_noe.tbl -i 7JU9_mc_cst.sc -n 500000 -s 7JU9_mc.sc -m m -e 0.025 -t 10.0 --outfile 7JU9_mc_ensemble.txt
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 7JVF_noe.tbl -i 7JVF_mc_cst.sc -n 500000 -s 7JVF_mc.sc -m m -e 0.025 -t 17.0 --outfile 7JVF_mc_ensemble.txt
```

Filter by atropisomers:
```
~/anaconda3/bin/python neg_ens/atripisomer_parse.py -p 1AJ1_mc_top.txt -d 17_CA-19_CA-15_CA-9_CA 9_CA-11_CA-7_CA-17_CA -r -f 1AJ1_ensemble.txt -o 1AJ1_mc_ensemble_filt.csv
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 1AJ1_noe.tbl -i 1AJ1_mc_cst.sc -n 500000 -s 1AJ1_mc.sc -m m -e 0.025 -t 10.0 -f 1AJ1_mc_ensemble_filt.csv --outfile 1AJ1_mc_ensemble.txt

~/anaconda3/bin/python neg_ens/atripisomer_parse.py -p 6PQG_mc_top.txt -d 16_CA-17_CA-14_CA-11_CA 15_CA-17_CA-13_CA-10_CA -r -f 6PQG_ensemble.txt -o 6PQG_mc_ensemble_filt.csv
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 6PQG_noe.tbl -i 6PQG_mc_cst.sc -n 500000 -s 6PQG_mc.sc -m m -e 0.025 -t 7.0 -f 6PQG_mc_ensemble_filt.csv --outfile 6PQG_mc_ensemble.txt
```

RMSD BB Calculation:
`echo "1AJ1 6PQG 6VHJ 6VLJ 7JU9 7JVF" | sed 's/ /\n&/g' | xargs -I % bash rmsd_bb_rosetta_mc.sh %`

Pymol session generation example:
```
item=1AJ1; pymol ~/Documents/EGFR_antibodies/RiPP_design/${item}.pdb $(cat ${item}_mc_ensemble.txt | grep -v sc | tr "\n" " ");
run prep_pymol.py
```

NMR Relax analysis
```
ls ????_free_nmr_relax.sc | xargs -I % echo % | sed 's/_free_nmr_relax.sc//' | while read -r base; do ~/anaconda3/bin/python ../gen_plots.py -i ${base}_free_nmr_relax.sc -c 10; done
ls 1AJ1_free_nmr_relax.sc 6PQG_free_nmr_relax.sc | xargs -I % echo % | sed 's/_free_nmr_relax.sc//' | while read -r base; do ~/anaconda3/bin/python ../gen_plots.py -i ${base}_free_nmr_relax.sc -c 30; done
ls ????_free_nmr_relax.sc | xargs -I % echo % | sed 's/_free_nmr_relax.sc//' | while read -r base; do bash plot_global_funnels.sh ${base}_free_nmr_relax_rmsd_bb.dat; done
```

