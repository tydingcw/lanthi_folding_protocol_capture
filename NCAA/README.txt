# NCAA Parameter Generation and Structure Prediction

This directory contains files and scripts related to the generation of parameters for non-canonical amino acids (NCAA) and their use in structure prediction.

## Directory Contents

- `fakerotlib_param_gen/`: Contains files generated and input into `fakerotlib` for parameter generation.
- `final_params/`: Stores the NCAA parameter files used for structure prediction, which are adopted from the outputs of `fakerotlib`.
- `gauss_opt/`: Scripts for geometry optimization, potential energy surface scanning, and generation of rama files for Rosetta.
- `sdf_files/`: Has the starting SDF files for the process.

## Scripts

- `make_3D.sh`: A script to convert 2D SDF files to 3D geometry.
- `min_bcl_conf.py`: Retrieves the best scoring conformer from BCL (Biochemical Library) sdf output.
- `run_molfile_parent.sh`: An example script demonstrating the use of `molfile_to_params_polymer.py` for parameter generation.


