#Clay Tydings

import argparse
import pandas as pd

def gen_ensemble_plot(score_file, ensemble_file):
    # Which structures do we care about (ensemble)
    ens_structs = []
    with open(ensemble_file, 'r') as in_file:
        for raw_line in in_file:
            line = raw_line.strip()
            if 'score' not in line:
                ens_structs.append(line.split('.')[0])

    # initialize the columns
    with open(score_file, 'r') as in_file:
        first_line = in_file.readline()
        first_line = in_file.readline()
        cols = first_line.split()

    out_str = ''

    for name in ens_structs:
        # for each file in ensemble, get the constraints
        df = pd.read_csv(score_file, delim_whitespace=True, skiprows=1)[cols]
        df = df.loc[df['description'] == name]
        #print(min(df['total_score'].tolist()), name) #TODO fix later
        temp_min = min(df['total_score'].tolist())
        temp_name = name.split('/')[-1]
        out_str += f'{temp_name}.pdb.gz {temp_min} \n'
    return out_str

parser = argparse.ArgumentParser(description = 'take in a tbl file and create rosetta nmr constraint file')
parser.add_argument('-e', '--ensemble', help = 'input ensemble list', required = True)
parser.add_argument('-s', '--scorefile', help = 'input scorefile', required = True)
parser.add_argument('-o', '--output', help = 'output name', required = True)
argument=parser.parse_args()

out_str = gen_ensemble_plot(argument.scorefile, argument.ensemble)

with open (argument.output, 'w') as outfile:
    outfile.write(out_str)