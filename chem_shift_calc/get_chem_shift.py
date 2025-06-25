import os
import argparse

def main():
    # for nmr shift need to take in bmr file
    parser = argparse.ArgumentParser(description='Read in chemical shifts from nmr experiment')

    # Add arguments
    parser.add_argument('-b', '--base', required=True,
                        help='The base of the output file.')
    parser.add_argument('-l', '--pdb_list', required=True,
                        help='The list of pdbs to iterate through.')
    parser.add_argument('-d', '--pdb_dir', required=True,
                        help='The directory with the pdbs to predict chemical shifts for.')

    args = parser.parse_args()
    pdb_dir = args.pdb_dir
    pdb_list = args.pdb_list
    #pdb_name = args.pdb_list.split('_')[0]
    pdb_name = args.base

    # Need a csv with chem_shift and description as columns
    csv_str = ''
    head = ''
    shift_str = ''
    index_count = 0
    first_read = True

    shift_dict = {}
    reading = False
    #print(f'{pdb_name}_shift.csv')
    with open(f'{pdb_name}_shift.csv', 'w') as shift_file:
        with open(pdb_list, 'r') as in_file:
            for raw_line in in_file:
                line = raw_line.strip()
                #split = line.split()
                #make the sparta predictions
                #print(f'gzip -c {pdb_dir}/{line} > temp.pdb ')
                os.system(f'gunzip -c {pdb_dir}/{line} > temp.pdb ')
                os.system('~/SPARTA+/SPARTA+/sparta+ -in temp.pdb')
                #extract data from pred.tab
                csv_str = ''
                shift_str = ''
                with open('pred.tab', 'r') as tab_file:
                    #for raw_tab in tab_file
                    #    tab = raw_tab.strip()
                    lines = tab_file.readlines()
                    vars_index = None
                    for i, data in reversed(list(enumerate(lines))): 
                        if 'VARS   RESID RESNAME ATOMNAME SS_SHIFT SHIFT RC_SHIFT HM_SHIFT EF_SHIFT SIGMA' in data:
                            vars_index = i
                            break
                    #print(lines[vars_index])
                    if vars_index != None:
                        for data_line in lines[vars_index + 3:]:
                            #0 is resi, 2 is atom, 4 is chem shift
                            #print(data_line)
                            split = data_line.split()
                            name = f'{split[0]}_{split[2]}'
                            shift = split[4]
                            #if csv_str == '':
                            csv_str += name + ','
                            shift_str += shift + ','
                        csv_str += 'description\n'
                        shift_str += line.replace('.pdb.gz', '').split('/')[1] + '\n'
                        if first_read:
                            head = csv_str
                            shift_file.write(csv_str)
                            shift_file.write(shift_str)
                        else:
                            if head != csv_str:
                                print('columns do match previous input')
                                exit(1)
                            else:
                                shift_file.write(shift_str)
                    else:
                        print(f'No VARS found in {line}')
                        exit(1)
                os.system(f'rm temp.pdb ')
                first_read = False
                #break
    #finished, report data 
    #print(shift_dict)

if __name__ == '__main__':
    main()
