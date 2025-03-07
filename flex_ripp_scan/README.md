# Flex RiPP Scan

This directory contains scripts and data for structure prediction and analysis of non-helical lanthipeptides.

The structure_prediction directory contains the scripts and data for repredicting lanthipeptide structures with Rosetta. 

## Scripts and Initial Setup

- Generate the initial XML and options for structure relaxation using `make_ripp_structpred_xml.py`:

  ```bash
  ~/anaconda3/bin/python ~/scripts/make_ripp_structpred_xml.py
  ```

- Use `sample_instruct.sh` to run Rosetta structure prediction scripts:

  ```bash
  ls options/6VHJ_* | xargs -I % -n1 -P 18 bash sample_instruct.sh %
  ls options/6V* | xargs -I % sbatch sample_instruct.slurm %
  ```

- Prepare options for a full run aiming for 200,000 structures using `fullrun_setup.py` and run `sample_instruct.sh` again:

  ```bash
  ~/anaconda3/bin/python fullrun_setup.py -i 6VLJ
  ls options_refined/6VHJ_*.options | xargs -I % -n1 -P 18 bash sample_instruct.sh %
  ```

- The `options_stored` directory was used to store successful options from iterative use of `make_ripp_structpred_xml.py`.

## Bluefin Job Array

Submit job arrays on Bluefin using the following commands:

```bash
sbatch --array=1-$(ls options/6VHJ_* | wc -l) sample_instruct_cmdarray.slurm options/6VHJ_ _rand_ind.options
# Repeat for other peptides as needed
```

## Post Processing

Analyze the results and generate plots:

```bash
~/anaconda3/bin/python ~/scripts/get_pert_rmsd_ddg.py -i 6VHJ_scan.sc -r rmsd_bb -s total_score -c 50.0
bash rmsd_native_pnear.sh 6VHJ_scan.dat
~/anaconda3/bin/python gen_plots.py -i 6VHJ_scan.sc
```

## Native Relax

Run native relax with the following command:

```bash
ls ????_f???.options ????_bb.options | xargs -n 1 -P 12 -I % bash -c 'bash sample_native.sh % `echo % | cut -c-4`'
```

## RMSD Calculations

Calculate RMSD for the peptides:

```bash
ls ????.cst | xargs -I % echo % | sed 's/.cst//' | while read -r base; do bash rmsd_bb.sh $base; done
```

## Peptide Modeling Information

Potential lanthipeptides to model and their characteristics:

- `1aj1`: 4 rings, 19 AA - DBB - acetonitrile/H2O
- `6pqg`: 4 rings, 19 AA - DBB - D2O/H2O
- `7ju9`: 3 rings, 20 AA - DBB - D2O/H2O
- `7jvf`: 2 rings, 21 AA - DBB - D2O/H2O
- `6vhj`: 2 rings, 16 AA - DBB - D2O/H2O - ProcA11
- `6vlj`: 2 rings, 19 AA - D2O/H2O - ProcA28

Note: `1AJ1` has a very old NOE constraint format. Replace `VAL QG qirh QQG` for residues 5 and 15 as needed.
```

