#import os
#import matplotlib.pyplot as plt
import argparse
#import numpy as np
from pymol import cmd
#from pymol import stored

def main():
    # We want to get a csv of all the relavent dihedral angles for each conformation of interest
    # need to iterate through the files, this can be done with the score file
    parser = argparse.ArgumentParser(description='Process iteration and channel arguments.')

    # Add arguments
    #parser.add_argument('-s', '--score', type=int, required=True,
    #                    help='The scorefile with descriptions.')
    parser.add_argument('-p', '--pdb_list', required=True,
                        help='The txt file with the pdb descriptions.')

    args = parser.parse_args()
    #score = args.score
    pdb_list = args.pdb_list
    pdb_name = args.pdb_list.split('_')[0]

    # Need a csv with phi/psi_resi and description as columns
    csv_str = ''
    resi_count = 0
    
    #print(f'{pdb_name}_dih.csv')
    with open(f'{pdb_name}_dih.csv', 'w') as dih_file:
        with open(pdb_list, 'r') as in_file:
            for raw_line in in_file:
                line = raw_line.strip()
                cmd.reinitialize()
                cmd.load(line)
                cmd.select('all')
                p = cmd.phi_psi('sele')
                #print(p)
                if resi_count == 0:
                    resi_count = len(p.keys())
                elif resi_count != len(p.keys()):
                    print('The number of residues in the peptide changed')
                    exit()
                if csv_str == '':
                    for i in range(2, resi_count+2):
                        csv_str += f'phi_{i},psi_{i},'
                    csv_str += 'description\n'
                    dih_file.write(csv_str)
                #add params here
                keys = list(p.keys())
                for i in range(0, resi_count):
                    #print(i)
                    phi = round(p[keys[i]][0], 2)
                    psi = round(p[keys[i]][1], 2)
                    #csv_str += f'{phi},{psi},'
                    dih_file.write(f'{phi},{psi},')
                #csv_str += line.replace('.pdb.gz', '') + '\n'
                dih_file.write(line.replace('.pdb.gz', '').split('/')[1] + '\n') #removes the extension and the directory
                #print(csv_str)
                #break

if __name__ == '__main__':
    main()
