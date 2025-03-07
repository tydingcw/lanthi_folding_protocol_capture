#import os
#import matplotlib.pyplot as plt
import argparse
#import numpy as np
#from pymol import cmd
#from pymol import stored

def main():
    # for nmr dih need to take in a talosn or mr(1aj1) file and csv of observed dihedrals
    parser = argparse.ArgumentParser(description='Read in dihedrals from nmr experiment')

    # Add arguments
    #parser.add_argument('-s', '--score', type=int, required=True,
    #                    help='The scorefile with descriptions.')
    parser.add_argument('-p', '--pred_dih', required=True,
                        help='The file with the backbone dihedral predictions.')

    args = parser.parse_args()
    #score = args.score
    pred_dih = args.pred_dih
    #pdb_name = args.pdb_list.split('_')[0]

    # Need a csv with phi/psi_resi and description as columns
    csv_str = ''
    resi_count = 0

    filetype = pred_dih.split('.')[-1]

    dih_dict = {}
    #print(f'{pdb_name}_dih.csv')
    #with open(f'{pdb_name}_dih.csv', 'w') as dih_file:
    if filetype == 'pred.tab': 
        print('pred.tab not yet supported')
        #with open(pred_dih, 'r') as in_file:
        #    for raw_line in in_file:
        #        line = raw_line.strip()
    elif filetype == 'mr': #1aj1
        reading = False
        with open(pred_dih, 'r') as in_file:
            for raw_line in in_file:
                line = raw_line.strip()
                if reading and line[0] == '!':
                    reading = False
                elif not reading:
                    if line == '#NMR_dihedral':
                        reading = True
                else: #reading is true and line is not !
                    split = line.split()
                    #print(split)
                    resi = split[1].split(':')[1].split('_')[1]
                    #resi = split[1]
                    #print(resi)
                    #check if it is a phi C-N-CA-C dihedral
                    if split[0][-1] == 'C' and split[1][-1] == 'N' and split[2][-2:] == 'CA' and split[3][-1] == 'C':
                        phi1 = float(split[4])
                        phi2 = float(split[5])
                        dih_dict[f'phi_{resi}'] = (min(phi1, phi2), max(phi1, phi2))
        #finished, report data 
        print(dih_dict)

if __name__ == '__main__':
    main()
