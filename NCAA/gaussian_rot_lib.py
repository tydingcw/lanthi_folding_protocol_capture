#code to make the rama tables for DHA
#Resname        phi     psi     probability     -log(prob)
#/home/tydingcw/Rosetta/main/database/scoring/score_functions/rama/alpha_ncaa/AIB_general.rama

#need to iterate over 0-35 and the r1, 
#0-35 represents the phi
#the psi angle is scanned
#double check that the dihedral is what is expected and that we have all 36
#make a dictionary of the phi, psi angles and the energies
#make a note of all combos that failed, assume highest observed energy
#calculate all probilities based on the lowest energy prob and normalize to 1
#then create the rama table based on AIB_general.rama

import numpy as np
import argparse
import scipy

parser = argparse.ArgumentParser(description = 'construct rama file from gaussian logs')
parser.add_argument('-l', '--log', help = 'example log file ending in _0.log', required = True)
argument = parser.parse_args()
log = argument.log
base = log.replace('_0.log', '_')
if base == log:
    print('bad log name')
name = base.split('_')[0]

rama_dict = {}
eng_vals = []
dih_list = [x*10 for x in range(36)]
failed_dih = []

for i in range(36):
    cnst_flag = False
    with open(f'{base}{i}.log', 'r') as in_file:
        curr_eng = 0
        phi = i*10
        psi = 0
        rama_dict[phi] = {}
        for line in in_file:
            if 'SCF Done:' in line:
                split = line.split()
                curr_eng = float(split[4])
            if 'Optimized Parameters' in line:
                rama_dict[phi][psi] = curr_eng
                eng_vals.append(curr_eng)
                #dih_vals.append(dih)
                #print('added '+str(phi)+' '+str(psi))
                psi += 10
            elif 'Error imposing constraints' in line or 'FormBX had a problem' in line:
                cnst_flag = True
                #print(f'{line} for {base}{i}.log')
    if cnst_flag:
        with open(f'{base}{i}_r1.log', 'r') as in_file:
            curr_eng = 0
            phi = i*10
            psi = 350
            for line in in_file:
                if 'SCF Done:' in line:
                    split = line.split()
                    curr_eng = float(split[4])
                if 'Optimized Parameters' in line:
                    if psi in rama_dict[phi].keys():
                        print(str(psi)+' duplicated')
                    rama_dict[phi][psi] = curr_eng
                    eng_vals.append(curr_eng)
                    #dih_vals.append(dih)
                    #print('added '+str(phi)+' '+str(psi))
                    psi -= 10
    print('scan '+str(phi)+' has ' + str(len(rama_dict[phi])) + ' psi angles')
    #print(phi, rama_dict[phi].keys())
    #print(dih_list)
    
    if len(rama_dict[phi]) < 36:
        failed_psi = [x for x in dih_list if x not in rama_dict[phi].keys()]
        #print(len(rama_dict[phi]), phi, len(failed_psi), failed_psi)
        for val in failed_psi:
            rama_dict[phi][val] = None
            failed_dih.append((phi, val))
    
min_eng = min(eng_vals)
max_eng = max(eng_vals)

#construct a scipy.interpolate object for missing values
phi_list = []
psi_list = []
eng_list = []
for phi_key in rama_dict.keys():
    for psi_key in rama_dict[phi_key].keys():
        if rama_dict[phi_key][psi_key] != None:
            phi_list.append(phi_key)
            psi_list.append(psi_key)
            eng_list.append(rama_dict[phi_key][psi_key])
interpolate = scipy.interpolate.SmoothBivariateSpline(phi_list, psi_list, eng_list)

#dealing with gaussian failure cases
for phi, psi in failed_dih:
    #rama_dict[phi][psi] = max_eng
    rama_dict[phi][psi] = interpolate(phi, psi)[0][0]

#units of energy are in hartrees (A.U.)
#pi/pj = e(Ej-Ei)/kT
#assume biological T is 37 Celsius, 310.15 kelvin
#T=310.15
T=300.0
#1 hartree is 4.3597E-18 J
hart_conv = 4.3597 * 10 ** -18
#k is 1.3806E-23 J/K
k = 1.3806 * 10 ** -23

#print(rama_dict)
rama_prob = {}
for i in range(36):
    rama_prob[i*10] = {}
    for j in range(36):
        if i*10 not in rama_dict.keys() or j*10 not in rama_dict[i*10].keys():
            print(i*10, j*10)
        else:
            rama_prob[i*10][j*10] = rama_dict[i*10][j*10] * hart_conv

#the sum of all prob variables will give the value to normalize by
tot_prob = 0
min_eng = min_eng * hart_conv
for i in range(36):
    for j in range(36):
        curr_eng = rama_prob[i*10][j*10]
        #print(min_eng-curr_eng)
        #print(min_eng, curr_eng)
        prob = np.exp((min_eng-curr_eng)/(k*T)) #prob is likihood of curr state compared to most likely state (min_eng)
        if min_eng == curr_eng:
            print(i, j)
        tot_prob += prob
        rama_prob[i*10][j*10] = prob
        
best_prob = (0, 0, 0)
#normalize probability
for i in range(36):
    for j in range(36):
        rama_prob[i*10][j*10] = rama_prob[i*10][j*10] / tot_prob
        if rama_prob[i*10][j*10] > best_prob[0]:
            best_prob = (rama_prob[i*10][j*10], i*10, j*10)
            
print(best_prob)

#write file (-180 to 170 degrees)
#Resname        phi     psi     probability     -log(prob)
rosetta_dih = [x*10-180 for x in range(36)]
with open(f'{name}.rama', 'w') as out_file:
    out_file.write('#Resname phi psi prob -log(prob)\n')
    for i in rosetta_dih:
        for j in rosetta_dih:
            if i < 0:
                i_new = 360 + i
            else:
                i_new = i
            if j < 0:
                j_new = 360 + j
            else:
                j_new = j
            prob = rama_prob[i_new][j_new]
            out_file.write(f'{name} {i} {j} {prob} {-np.log(prob)}\n')
