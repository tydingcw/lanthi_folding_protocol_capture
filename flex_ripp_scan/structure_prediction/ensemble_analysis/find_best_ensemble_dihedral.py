##!/bin/env python2.7
#Adapted from "Integrating Solid-State NMR and Computational Modeling to Investigate the Structure and Dynamics of Membrane-Associated Ghrelin"
#ghrelin/scripts/python_scripts

#from optparse import OptionParser
import argparse
import subprocess
import math
import glob
import numpy as np
import random
import sys
import pandas as pd
import os

# functions
# data set parsing, take quantile 0.2 lowest energy and get their measured distances
def parse_file(filename, scorefile, energy=0.2, energy_list = [], pred_dih=None, obs_dih=None):
	structs = [] # these are the structs passing the filter
	base = ''
	data = {}
	#currently, this pairs (res,atom) with (expCS,predCS)
	#we just need to look at all of the nmr constraints ex. 8CWX_cst.sc score, measure

	# probably want to limit to low energy structures, top 20% or so
	df = pd.read_csv(scorefile, delim_whitespace=True, skiprows=1)
	score_threshold = df['total_score'].quantile(energy)
	print(f'Energy threshold is: {score_threshold}')
	df = df[df['total_score'] <= score_threshold]
	if (energy_list == [] or pred_dih==None or obs_dih==None) : 
		for desc in df['description']: #TODO add filter for backbone dihedral
			#print('desc', desc)
			structs.append(desc.split('/')[-1])
			base = desc.split('/')[0]
		#print(len(structs))
	else:
		#we are refining this with dihedral constraints
		for desc in energy_list:
			dih_score = score_sing_dih(pred_dih, obs_dih[desc])
			if dih_score < 1000:
				structs.append(desc)
		print(f'passing dih filter: {len(structs)}')

	#probably want to compare the measured difference - will use nmr formula
	with open(filename) as in_file:
		first_line = in_file.readline()
		first_line = in_file.readline()
		cols=first_line.split()
	df = pd.read_csv(filename, delim_whitespace=True, skiprows=1)[cols]

	for desc in df['description']:
		#structs.append(desc)
		if desc in structs:
			#print(desc)
			df_temp = df.loc[df['description'] == desc ]
			data[desc] = {}
			for col in cols:
				if 'measure' in col and 'CstFile' in col:
					#print(col.split('_')[1], df_temp[col])
					num = int(col.split('_')[1])
					measure = float(df_temp[col])
					data[desc][num] = measure

		#for index,line in enumerate(filehandle):
		#	if index == 0:
		#		continue
		#	if line.startswith("#yes"):
		#		return (filename,data)
		#	line = line.split()
		#	res = int(line[0])
		#	atom = line[2]
		#	expCS = float(line[3])
		#	predCS = float(line[7])
		#	data[(res,atom)] = (expCS,predCS)
		#	print data
		#	if found_header:
		#		line = line.split()
	return data, base

def parse_negative(negative, structs):
	#structs = []
	data = {}

	#probably want to compare the measured difference - will use nmr formula
	with open(negative) as in_file:
		first_line = in_file.readline()
		first_line = in_file.readline()
		cols=first_line.split()
	df = pd.read_csv(negative, delim_whitespace=True, skiprows=1)[cols]

	for desc in df['description']:
		#structs.append(desc)
		if desc in structs:
			#print(desc)
			df_temp = df.loc[df['description'] == desc ]
			data[desc] = {}
			for col in cols:
				if 'measure' in col and 'CstFile' in col:
					#print(col.split('_')[1], df_temp[col])
					num = int(col.split('_')[1])
					measure = float(df_temp[col])
					data[desc][num] = measure

	return data

def parse_obs_dih(obs_dih, structs):
	data = {}

	with open(obs_dih) as in_file:
		first_line = in_file.readline().strip()
		cols=first_line.split(',')
	#print(cols)
	df = pd.read_csv(obs_dih, delim_whitespace=False)[cols]
	#df = pd.read_csv(obs_dih, delim_whitespace=False)
	#print(df)

	for desc in df['description']:
		if desc in structs:
			#print(desc)
			df_temp = df.loc[df['description'] == desc ]
			data[desc] = {}
			for col in cols:
				if col != 'description':
					measure = float(df_temp[col])
					data[desc][col] = measure

	return data
def parse_pred_dih(pred_dih):
    filetype = pred_dih.split('.')[-1]

    dih_dict = {}
    #print(f'{pdb_name}_dih.csv')
    #with open(f'{pdb_name}_dih.csv', 'w') as dih_file:
    if filetype == 'pred.tab': 
        print('pred.tab not yet supported')
        exit(1)
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
        #print(dih_dict)
    else:
        print(f'{pred_dih} is not supported')
        exit(1)
    return dih_dict

#gets the center distances from the nmr
def parse_reference(filename):
	ref_dict = {}
	aa_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
           'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
           'DBB', 'DBU', 'DHA', 'DBS', 'DBR', 'DAL', 'ABA'] #'BDS', 'BLR'
	with open (filename, 'r') as infile:
		count = 0
		for raw_line in infile:
			line = raw_line.strip()
			#if argument.fix:
			#	line = fix_naming(line)
			split = line.split()
			hash_split = line.split('#')[0].split()
			if len(line) == 0 or line[0] == '#':
				pass
			elif line[0:6] == 'assign': #this is for XPLOR type files
				#res1 = str(int(split[2]) + start)
				#res2 = str(int(split[7]) + start)
				#atom1 = split[5][:-1]
				#atom2 = split[10][:-1]
				center = float(split[11])
				#lower = float(split[12])
				#upper = float(split[13])
				#ambig = False
				#print(line)
				#atom1, atom2 = fix_conn_naming(res1, res2, atom1, atom2)
				if center > 2.3 :
					count += 1
					ref_dict[count] = center
			elif len(split) == 64 and split[8] in aa_list and split[18] in aa_list: #2KTO
				#res1 = str(int(split[7]) + start)
				#res2 = str(int(split[17]) + start)
				#atom1 = split[51]
				#atom2 = split[59]
				center = float(split[28])
				#lower = center - 1.8
				#upper = 0.0
				#lets fix any connection naming
				#atom1, atom2 = fix_conn_naming(res1, res2, atom1, atom2)
				if center > 2.3 :
					#temp_str = make_constraint(res1, res2, atom1, atom2, center, lower, upper)
					count += 1
					ref_dict[count] = center
					#out_str += temp_str.replace('count', str(count))
			elif len(hash_split) == 7 and split[1] in aa_list and split[4] in aa_list:  # 6PQG
				center = float(split[6])
				if center > 2.3:
					count += 1
					ref_dict[count] = center
			else:
				print(f'{split[8]} or {split[18]} not an AA')
	return ref_dict

# monte carlo moves
def mc(ncycles,dataset,reference,mode,temp,negative=None, neg_cutoff=4.0, pred_dih=None, obs_dih=None):
	current_dataset = [] #dataset
	for key in dataset:
		current_dataset.append((key,dataset[key])) #create list with tag,data
	last_dataset = current_dataset #dataset
	min_ensemble_size = args.min_ensemble_size
	max_ensemble_size = args.max_ensemble_size
	current_ensemble = []
	last_ensemble = current_ensemble
	best_ensemble = current_ensemble
	last_score = 1000000.0
	best_score = 1000000.0
	accept_count = 0.0
	reject_count = 1.0
	accept_ratio = 0.0

	# add models from the pool to get a minimum ensemble
	for i in range(ncycles):
		while len(current_ensemble) < min_ensemble_size:
			add(current_dataset,current_ensemble)
		with open('best_score.txt', 'a') as dih_file:
			dih_file.write(f'{best_score}\n')
		# 1. make some move
		if len(current_ensemble) > min_ensemble_size and len(current_ensemble) < max_ensemble_size:
			# functions are objects in python
			moves = [add,remove,swap]
			move = random.choice(moves)
			move(current_dataset,current_ensemble)
		elif len(current_ensemble) == min_ensemble_size:
			moves = [add,swap]
			move = random.choice(moves)
			move(current_dataset,current_ensemble)
		elif len(current_ensemble) == max_ensemble_size:
			moves = [remove,swap]
			move = random.choice(moves)
			move(current_dataset,current_ensemble)
		# 2. score move
		if negative != None:
			current_score = score_neg(current_ensemble,reference, negative, neg_cutoff)
		elif pred_dih != None and obs_dih != None:
			current_score = score_dih(current_ensemble,reference, pred_dih, obs_dih)
		else:
			current_score = score(current_ensemble,reference)
		# 3. accept or reject based on metropolis
		if mode == 'o': #keeping this as teh default, original script mode
			if current_score <= last_score:
				# accept move
				last_score = current_score
				last_ensemble = current_ensemble
				last_dataset = current_dataset
				accept_count += 1.0
			else:
				# reject move
				current_score = last_score
				current_ensemble = last_ensemble
				current_dataset = last_dataset
				reject_count += 1.0

		elif mode == 'm':
			if current_score <= last_score:
				# accept move
				last_score = current_score
				last_ensemble = current_ensemble
				last_dataset = current_dataset
				accept_count += 1.0
			else:
				boltzmann = (last_score - current_score) / temp
				probability = math.exp(boltzmann)
				rand = random.uniform(0.0,1.0)
				#print(f'last_score: {last_score}, current_score: {current_score}, boltz {boltzmann}, prob: {probability}, temp: {temp}, rand: {rand}')
				#if probability >= 1.0:
				#	print('exiting because probability is 1 or greater')
				#	print(f'last_score: {last_score}, current_score: {current_score}, boltz {boltzmann}, prob: {probability}')
				#	exit(1)
				if probability < 1.0 and probability < rand: #prob closer to one means lower energy gap, should be higher acceptance
					# reject move
					#print('reject based on metropolis')
					current_score = last_score
					current_ensemble = last_ensemble
					current_dataset = last_dataset
					reject_count += 1.0
				else:
					# accept move
					#print('accept based on metropolis')
					last_score = current_score
					last_ensemble = current_ensemble
					last_dataset = current_dataset
					accept_count += 1.0
		else:
			print('mode not supported')
			exit(1)
		if last_score <= best_score:
			best_score = last_score
			best_ensemble = last_ensemble
		# print accept ratio (for optimization)
		accept_ratio = accept_count / (accept_count + reject_count)
		# print("cycle:  " + str(i) + "\tbest_score:  " + str(best_score) + "\tacceptance_ratio:  " + str(accept_ratio))
	return (best_ensemble,best_score,accept_ratio)

# objective function (average RMSD of the ensemble decreases)
def score(current_set,reference): # current_set is a list of maps (filename, map)
	# get experimental values
	first_record = current_set[0][1] # map value is (expCS, predCS)
	#exp_map replaced by reference
	#exp_map = {}
	#for key in first_record:
	#	exp_map[key] = first_record[key][0]
	#print str(exp_map)+"\n"
	#print "len(exp_map):  " + str(len(exp_map)) + "\n"
	# Now get the predicted values for each given key (experimental CS)
	pred_map = {}
	diffs = []
	for key in reference: # this chould be 1,2,... for all nmr_cst numbers
		pred_values = []
		for filename,data in current_set:
			pred_values.append(data[key])
		#print('debug')
		#print(pred_values)
		#pred_map[key] = np.mean(pred_values)
		#print(pred_map[key])
		#print(np.array(pred_values) ** -6)
		#print((np.sum(np.array(pred_values) ** -6)))
		#print((np.sum(np.array(pred_values) ** -6)/ len(pred_values))**(-1/6))
		pred_map[key] = (np.sum(np.array(pred_values) ** -6)/ len(pred_values))**(-1/6)
		score = abs(reference[key] - pred_map[key]) ** 2
		diffs.append(score)
#	print str(pred_map) + "\n"
#	print "len(pred_map):  " + str(len(pred_map)) + "\n"
	#want (1/N Sum dist^-6)^-1/6
	#for key in pred_map:
		#diffs.append((pred_map[key] - reference[key])**2)
		#diffs.append(pred_map[key] - reference[key])
		#diffs.append((pred_map[key])**-6)
	#rmsd = numpy.sqrt(numpy.sum(diffs)/len(diffs))
	sum = np.sum(diffs)
	#rmsd = numpy.sqrt(sum)
	#rmsd = abs(reference[key] - sum**(-1/6))
	#print(diffs)
	#print(rmsd)
	return sum

def score_neg(current_set, reference, negative, cutoff=4.0):  # current_set is a list of maps (filename, map)
	#pred_map = {}
	#diffs = []
	sum = 0
	neg_map = {}
	#for key in reference:  # this chould be 1,2,... for all nmr_cst numbers
	#	pred_values = []
	#	for filename, data in current_set:
	#		pred_values.append(data[key])
	#	pred_map[key] = (np.sum(np.array(pred_values) ** -6)/ len(pred_values))**(-1/6)
	#	score = abs(reference[key] - pred_map[key]) ** 2
	#	#diffs.append(score)
	#	sum += score
	#	#pred_map[key] = np.mean(pred_values)
	sum = score(current_set, reference)
	#print(f'score {sum}')
	#print(f'non neg score {sum} {score(current_set, reference)}')
	#diffs = []
	# want (1/N Sum dist^-6)^-1/6
	#for key in pred_map:
	#	diffs.append((pred_map[key]) ** -6)
	#sum = np.sum(diffs) #/ len(diffs)
	#rmsd = abs(reference[key] - sum ** (-1 / 6))
	#iterate through negative
	#score = 0
	first_struct = list(negative.keys())[0]
	#print(list(negative[first_struct])[0:10])
	#violations = 0
	for key in list(negative[first_struct]): #all constraints:
		pred_values = []
		for filename, data in current_set:
			pred_values.append(negative[filename][key])
		# want (1/N Sum dist^-6)^-1/6
		neg_map[key] = (np.sum(np.array(pred_values) ** -6)/ len(pred_values))**(-1/6)
		#sum = np.sum(neg_map[key])/len(neg_map[key])
		#score = min(5.0 - sum ** (-1/6), 0) ** 2
		diff = min(neg_map[key] - cutoff, 0) ** 2
		#if diff > 0:
		#	violations +=1
		sum += diff
	#print(f'neg score {sum}')
	#print(f'{len(list(negative[first_struct]))} {violations}')
	return sum

def score_dih(current_set, reference, pred_dih, obs_dih):  # current_set is a list of maps (filename, map)
	sum = 0
	neg_map = {}
	sum = score(current_set, reference)

	dih_score = 0
	first_struct = list(obs_dih.keys())[0]
	for key in [x for x in list(pred_dih.keys()) if x in obs_dih[first_struct].keys()]: #all dihedrals except first and last
		#pred_values = []
		for filename, data in current_set:
			#get the dihedral value
			dih_val = obs_dih[filename][key]
			deviation = 0
			if dih_val >= pred_dih[key][0] and dih_val <= pred_dih[key][1]:
				pass #deviation is 0
			else:
				#need to see how far outside of the range we are
				#https://gamedev.stackexchange.com/questions/4467/comparing-angles-and-working-out-the-difference
				dev1 = 180 - abs(abs(dih_val - pred_dih[key][0]) - 180); 
				dev2 = 180 - abs(abs(dih_val - pred_dih[key][1]) - 180); 
				deviation = min(dev1, dev2)
			diff = deviation ** 2
			dih_score += diff
			#pred_values.append(negative[filename][key])
		#neg_map[key] = (np.sum(np.array(pred_values) ** -6)/ len(pred_values))**(-1/6)
		#diff = min(neg_map[key] - cutoff, 0) ** 2
		#sum += diff
	with open('dih_file.txt', 'a') as dih_file:
		dih_file.write(f'{dih_score}\n')
	#sum + dih_score
	return sum

def score_sing_dih(pred_dih, obs_dih):
	dih_score = 0
	for key in [x for x in list(pred_dih.keys()) if x in obs_dih.keys()]: #all dihedrals except first and last
		dih_val = obs_dih[key]
		deviation = 0
		if dih_val >= pred_dih[key][0] and dih_val <= pred_dih[key][1]:
			pass #deviation is 0
		else:
			dev1 = 180 - abs(abs(dih_val - pred_dih[key][0]) - 180)
			dev2 = 180 - abs(abs(dih_val - pred_dih[key][1]) - 180)
			deviation = min(dev1, dev2)
		diff = deviation ** 2
		dih_score += diff
	return dih_score

#	data_for_scoring = [x[1] for x in all_data] # for x in all_data, append x[1] to data_for_scoring
#	return numpy.mean(data_for_scoring)


# add to ensemble
def add(pool,test):
	item_to_add = random.sample(pool,1)[0]
	test.append(item_to_add)
	pool.remove(item_to_add)

# remove from ensemble
def remove(pool,test):
	item_to_remove = random.sample(test,1)[0]
	pool.append(item_to_remove)
	test.remove(item_to_remove)

# swap from ensemble
def swap(pool,test):
	test_random_item = random.sample(test,1)[0]
	pool_random_item = random.sample(pool,1)[0]
	pool.append(test_random_item)
	test.append(pool_random_item)
	pool.remove(pool_random_item)
	test.remove(test_random_item)

# options stuff
usage="%prog [options] <input_directory>"
#parser=OptionParser(usage)
parser = argparse.ArgumentParser(description = 'take in tbl file and scored cnst file, find best ensemble for cst')
parser.add_argument('-i', '--input', help = 'input score file with constraints', required = True)
parser.add_argument('-r', '--reference', help = 'reference tbl file for optimal distance', required = True)
parser.add_argument('-a', '--neg_ref', help = 'reference file for negative nmr data to prevent interactions', default = '')
parser.add_argument('-c',"--neg_cutoff",help="cutoff to use for negative noe",type=float,default=5.0)
parser.add_argument('-s', '--scorefile', help = 'scorefile to filter high energy struct based on percentile', required = True)
parser.add_argument('-t',"--temp",help="temperature to use for MMC cycles",type=float,default=0.01)
parser.add_argument('-n', '--ncycles',dest="ncycles",help="number of MMC cycles",type=int,default=10000)
#potential to add simulated annealing(want temp to start low enough, lower than 1(log(1+i)))
parser.add_argument('-m', '--mode',help="which type of simulation to perform o-optimise, m-metropolis",default="m")
parser.add_argument("--min_ensemble_size",dest="min_ensemble_size",help="mininum ensemble size",type=int,default=10)
parser.add_argument("--max_ensemble_size",dest="max_ensemble_size",help="maximum ensemble size",type=int,default=30)
parser.add_argument("--outfile",dest="outfile",help="outfile",default="")
parser.add_argument('-e', "--energy",help="percentile cutoff for energy", type=float, default=0.2)
# for nmr dih need to take in a talosn or mr(1aj1) file and csv of observed dihedrals
parser.add_argument("--pred_dih",dest="pred_dih",help="file with predicted dih from NMR CS, only support 1aj1.mr type files right now",default="")
parser.add_argument("--obs_dih",dest="obs_dih",help="file with observed dih for the top quantile of structures",default="")
args = parser.parse_args()
#(opts,args) = parser.parse_args()

# set up variables
ncycles = args.ncycles
#temp = float(opts.temp)
#dataset = []
if args.mode not in ['o', 'm']:
	print(f'{args.mode} is not a supported mode')
	exit(1)

# function calls
# parse the files
#print "parsing files..."
#input_directory = args[0]
#this would create a dataset of reference values and calcuated values
#for file in glob.glob("%s/*.cs.out" %input_directory):
#dataset.append(parse_file(args.input))
dataset, base = parse_file(args.input, args.scorefile, energy=args.energy)
#print(len(dataset.keys()))
negative = None
if args.neg_ref != '':
	negative = parse_negative(args.neg_ref, dataset.keys())

if (args.pred_dih != '' and args.obs_dih == '') or (args.pred_dih == '' and args.obs_dih != ''):
    print('both pred_dih and obs_dih should be defined or undefined')
    exit(1)

pred_dih_vals = None
if args.pred_dih != '':
    pred_dih_vals = parse_pred_dih(args.pred_dih)

obs_dih_vals = None
if args.obs_dih != '':
    obs_dih_vals = parse_obs_dih(args.obs_dih, dataset.keys())
#print "...done."
#print(len(negative.keys()))
#print(len([x for x in dataset.keys() if x not in negative.keys()]))
reference = parse_reference(args.reference)

# run the MC cycles
print(f'starting {ncycles} steps')
if args.neg_ref != '':
	mc_result = (mc(ncycles,dataset,reference,args.mode,args.temp,negative,args.neg_cutoff))
elif args.pred_dih != '' and args.obs_dih != '':
	#TODO refilter dataset based on dihedrals
	dataset, temp = parse_file(args.input, args.scorefile, energy=args.energy, energy_list=dataset.keys(), pred_dih=pred_dih_vals, obs_dih=obs_dih_vals)
	mc_result = (mc(ncycles,dataset,reference,args.mode,args.temp, pred_dih=pred_dih_vals, obs_dih=obs_dih_vals))
else:
	mc_result = (mc(ncycles,dataset,reference,args.mode,args.temp))
best_ensemble = mc_result[0]
best_score = mc_result[1]
accept_ratio = mc_result[2]
print(f'best score: {best_score}')
cst_num = len(reference.keys())
print(f'num cst: {cst_num}')
print(f'avg cst score: {best_score/cst_num}')
print(f'accept ratio: {accept_ratio}')

# print ensemble of models
for item in best_ensemble:
	print(f'{base}/{item[0]}.pdb.gz')
if args.outfile == '':
	if os.path.exists(f'{base}_ensemble.txt'):
		print(f'{base}_ensemble.txt already exits')
		for item in best_ensemble:
			print(f'{base}/{item[0]}.pdb.gz')
	else:
		with open(f'{base}_ensemble.txt','w') as outfile:
			outfile.write("score:  "+str(best_score)+"\n")
			for item in best_ensemble:
				tag =f'{base}/{item[0]}.pdb.gz'
				outfile.write(tag+"\n")
else:
	with open(args.outfile,'w') as outfile:
		outfile.write("score:  "+str(best_score)+"\n")
		for item in best_ensemble:
			tag =f'{base}/{item[0]}.pdb.gz'
			outfile.write(tag+"\n")
