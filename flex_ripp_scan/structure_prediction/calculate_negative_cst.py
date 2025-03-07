import argparse
import gzip
import os

def parse_pdb(pdb_file):
    pdb_data = {}
    if '.gz' in pdb_file:
        with gzip.open(pdb_file, 'rb') as file:
            for raw_line in file:
                line = raw_line.decode('utf-8')
                split = line.split()
                if (line.startswith('ATOM') or line.startswith('HETATM')) and split[-1] == 'H':
                    residue_position = int(line[22:26].strip())
                    atom_name = line[12:16].strip()
                    if residue_position not in pdb_data:
                        pdb_data[residue_position] = []
                    if atom_name not in pdb_data[residue_position]:
                        pdb_data[residue_position].append(atom_name)
    else:
        with open(pdb_file, 'r') as file:
            for line in file:
                split = line.split()
                if (line.startswith('ATOM') or line.startswith('HETATM')) and split[-1] == 'H':
                    residue_position = int(line[22:26].strip())
                    atom_name = line[12:16].strip()
                    if residue_position not in pdb_data:
                        pdb_data[residue_position] = []
                    if atom_name not in pdb_data[residue_position]:
                        pdb_data[residue_position].append(atom_name)
    return pdb_data

def parse_cst(cst_file):
    cst_combos = set()
    with open(cst_file, 'r') as file:
        for line in file:
            split = line.split()
            resi1 = int(split[2])
            resi2 = int(split[4])
            cst_combos.add((resi1, resi2))
            cst_combos.add((resi2, resi1))
    return cst_combos

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Create a new constraint list from a PDB file.')
    parser.add_argument('-p', '--pdb_file', type=str, help='Path to the PDB file', required=True)
    parser.add_argument('-c', '--cst_list', type=str, help='Path to the constraint list file', required=True)
    args = parser.parse_args()

    # Check if constraint list file exists
    out_name = args.cst_list.replace('.cst', '') + '_neg.cst'
    if os.path.exists(out_name):
        print(f"Constraint list file '{out_name}' already exists. No new file will be created.")
    else:
        # Generate the new constraint list from the PDB file
        neg_cst = ''
        resi_dict = parse_pdb(args.pdb_file)
        cst_set = parse_cst(args.cst_list)
        #print(cst_set)
        resi_list = list(resi_dict.keys())
        for resi in resi_list:
            for partner in [x for x in resi_list if x < resi]:
                if ((resi, partner) not in cst_set) and ((partner, resi) not in cst_set):
                    #print(resi, partner)
                    for resi_hyd in resi_dict[resi]:
                        for part_hyd in resi_dict[partner]:
                            neg_cst += f'AtomPair {resi_hyd} {resi} {part_hyd} {partner} HARMONIC 5.0 1.0 \n'
        with open(out_name, 'w') as out_file:
            out_file.write(neg_cst)
        print(f"New constraint list file '{out_name}' created.")

if __name__ == "__main__":
    main()
