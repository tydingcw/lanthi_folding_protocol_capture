import argparse
import pandas as pd
import os

# functions
# data set parsing, take quantile 0.2 lowest energy and get their measured distances
def parse_file(scorefile, energy=0.2):
	structs = ''
	base = ''
	# probably want to limit to low energy structures, top 20% or so
	df = pd.read_csv(scorefile, delim_whitespace=True, skiprows=1)
	score_threshold = df['total_score'].quantile(energy)
	print(f'Energy threshold is: {score_threshold}')
	df = df[df['total_score'] <= score_threshold]
	for desc in df['description']:
		#print('desc', desc)
		structs += desc + '.pdb.gz\n'

	return structs

# options stuff
usage="%prog [options] <input_directory>"
#parser=OptionParser(usage)
parser = argparse.ArgumentParser(description = 'take in rosetta scorefile and output the top percentile of structures')
parser.add_argument('-i', '--input', help = 'input score file with constraints', required = True)
parser.add_argument("--outfile",dest="outfile",help="outfile",default="")
parser.add_argument('-e', "--energy",help="percentile cutoff for energy", type=float, default=0.025)
args = parser.parse_args()


dataset = parse_file(args.input, energy=args.energy)

base = args.input.split('_')[0]

# print ensemble of models
if args.outfile == '':
	if os.path.exists(f'{base}_top.txt'):
		print(f'{base}_top.txt already exits')
	else:
		with open(f'{base}_top.txt','w') as outfile:
			outfile.write(dataset)
else:
	with open(args.outfile,'w') as outfile:
		outfile.write(dataset)
