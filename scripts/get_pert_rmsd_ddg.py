#import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description = 'take in a score file and output a csv with rmsd and total_score')
parser.add_argument('-i', '--in_file', help = 'input *.sc', required = True)
parser.add_argument('-r', '--rmsd', help = 'rmsd name', default='rmsd')
parser.add_argument('-s', '--score', help = 'score type', default='ddg_bind')
parser.add_argument('-c', '--cutoff', help = 'score cutoff', type=float, default=0.0)

argument=parser.parse_args()

#filename = sys.argv[1]

cols = None
with open(argument.in_file, 'r') as in_file:
    first_line = in_file.readline()
    first_line = in_file.readline()
    cols=first_line.split()

#print(cols)

df = pd.read_csv(argument.in_file, delim_whitespace=True, skiprows=1)[cols]
#df = df[['ddg_bind','rmsd']]
df = df[[argument.score,argument.rmsd]]
df = df[df[argument.score] < argument.cutoff]
#print(df.columns)
#for des in df.nsmallest(n=argument.number, columns=[argument.column])['description']:
#    print(des)

#out_name = argument.in_file.replace('.sc', f'_{argument.rmsd}.dat')
out_name = argument.in_file.replace('.sc', '.dat')
print(out_name)
df.to_csv(out_name, index=False, sep = ' ')
