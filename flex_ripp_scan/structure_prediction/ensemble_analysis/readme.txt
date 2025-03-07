find_best_ensemble.py should find the best ensemble using nmr distance rmsd

Need to get the correct distances to compare to:
~/scripts/make_nmr_cst.py shows how to parse

~/anaconda3/bin/python find_best_ensemble.py -r ../2KTO_noe.tbl -i ../mem_try_final/2KTO_cst.sc -n 10000 -s ../mem_try_final/2KTO_af_sample.sc

Generating Ensemble:
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 2KTO_noe.tbl -i 2KTO_cst.sc -n 500000 -s 2KTO_af_sample.sc -m m -t 0.01
#actual cmd:
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 2KTO_noe.tbl -i 2KTO_cst.sc -n 500000 -s 2KTO_af_sample.sc -m m -t 0.005
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 8CWX_noe.tbl -i 8CWX_cst.sc -n 500000 -s 8CWX_af_sample.sc -m m -t 0.005
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 6VGT_noe.tbl -i 6VGT_cst.sc -n 500000 -s 6VGT_af_sample.sc -m m -t 0.01
~/anaconda3/bin/python ensemble_analysis/find_best_ensemble.py -r 6VE9_noe.tbl -i 6VE9_cst.sc -n 500000 -s 6VE9_af_sample.sc -m m -t 0.02

cp ../../nmr_cst/????_noe.tbl ./

