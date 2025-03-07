# AlphaFold Predictions for Helical Lanthipeptides

This directory contains the AlphaFold prediction data for helical lanthipeptides and further analysis of these peptides.

## Directory Contents

- `2KTO/`, `6VE9/`, `6VGT/`, `8CWX/`: Each directory contains the AlphaFold prediction models for the respective peptide.
- `2KTO.fasta`, `6VE9.fasta`, `6VGT.fasta`, `8CWX.fasta`: FASTA files with the amino acid sequences used for the AlphaFold predictions.
- `run_alphafold_ACCRE_Turing.sb`: Slurm batch script to run AlphaFold predictions on the ACCRE Turing cluster.
- `top_pred/`: Directory where the top-ranked AlphaFold predictions with lanthionine rings and other non-canonical amino acids (NCAA) introduced are stored. The structures of these peptides are predicted with Rosetta, and additional analysis comparing Rosetta predicted structures and PDB structures is included.

## Running AlphaFold Predictions

To run AlphaFold predictions, ensure that all prerequisite software and environment settings are properly configured on your system, then submit the job to Slurm with the following command:

```bash
sbatch run_alphafold_ACCRE_Turing.sb fasta_file.fasta

setting early template date, use all 5 models

for file in 2KTO/ranked_?.pdb; do cp $file top_pred/2KTO_$( basename $file | sed 's/ranked_//' ); done
