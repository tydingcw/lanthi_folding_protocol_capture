import argparse

parser = argparse.ArgumentParser(description = 'take in a tbl file and create rosetta nmr constraint file')
parser.add_argument('-i', '--input', help = 'input tbl file', required = True)
parser.add_argument('-c', '--conn', help = 'file specifying lanthionine connectivity, needed to fix atom naming in some cases', default='')
parser.add_argument('-s', '--seq', help = 'input sequence file', required = True)
parser.add_argument('-l', '--loose', help = 'should nmr constraints be exact or loose', type=int, default = 1)
parser.add_argument("-f","--fix",action="store_true",  help="Fix naming for input. i.e. remove HB3 for CYS. Be careful with this. ", default=False) #superseeded by -name
parser.add_argument('-n', '--name', help="if AA naming is not HB1 and HB2, is HB2 and HB3, specify not-norm here", default='norm')
#parser.add_argument('-n', '--nmr', help = 'name of the nmr constraint type (xplor, dyana), default is', default = 'xplor')
argument=parser.parse_args()

name_dict = {'CYS': {'H': 'H', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB'}, 
             'HIS': {'HE1': 'HE1', 'HE2': 'HE2', 'H': 'H', 'HB1': '1HB', 'HA': 'HA', 'HB2': '2HB', 'HD2': 'HD2'}, 
             'GLY': {'H': 'H', 'HA2': '2HA', 'HA1': '1HA'}, 
             'DAL': {'H': 'H', 'HA': 'HA', 'HB3': '3HB', 'HB2': '2HB', 'HB1': '1HB'}, 
             'SER': {'H': 'H', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB', 'HG': 'HG'}, 
             'PRO': {'HD2': '2HD', 'HD1': '1HD', 'HG1': '1HG', 'HG2': '2HG', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB',}, 
             'ALA': {'H2': '2H', 'H3': '3H', 'H1': '1H', 'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB3': '3HB', 'HB2': '2HB'}, 
             'DHA': {'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB2': '2HB'}, 
             'DBU': {'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB2': '2HB', 'HG13': '3HG', 'HG12': '2HG', 'HG11': '1HG'}, 
             'DALA': {'H2': '2H', 'H3': '3H', 'H1': '1H', 'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB3': '3HB', 'HB2': '2HB'}, 
             'MET': {'HE1': '1HE', 'HE2': '2HE', 'HE3': '3HE', 'H': 'H', 'HG1': '1HG', 'HG2': '2HG', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB'}, 
             'GLU': {'H': 'H', 'HG1': '1HG', 'HG2': '2HG', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB'}, 
             'GLN': {'H': 'H', 'HG1': '1HG', 'HG2': '2HG', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB', 'HE22': '2HE2', 'HE21': '1HE2'}, 
             'TRP': {'HH2': 'HH2', 'HE1': 'HE1', 'HD1': 'HD1', 'HE3': 'HE3', 'H': 'H', 'CH2': 'CH2', 'HZ3': 'HZ3', 'HZ2': 'HZ2', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB'}, 
             'ASN': {'HD22': '2HD2', 'HD21': '1HD2', 'H': 'H', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB'}, 
             'ASP': {'H': 'H', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB'}, 
             'TYR': {'HE2': 'HE2', 'HD2': 'HD2', 'HD1': 'HD1', 'OH': 'OH', 'H': 'H', 'HE1': 'HE1', 'HH': 'HH', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB'}, 
             'PHE': {'HE2': 'HE2', 'HD2': 'HD2', 'HD1': 'HD1', 'HZ': 'HZ', 'H': 'H', 'HE1': 'HE1', 'HH': 'HH', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB'}, 
             'VAL': {'HG1': 'HG1', 'HG2': 'HG2', 'H': 'H', 'HA': 'HA', 'HB': 'HB', 'HG11': '1HG1', 'HG12': '2HG1', 'HG13': '3HG1', 'HG21': '1HG2', 'HG22': '2HG2', 'HG23': '3HG2'}, 
             'ARG': {'HH21': '1HH2', 'HH22': '2HH2', 'HH11': '1HH1', 'HE12': '2HH1', 'HE': 'HE', 'HD1': '1HD', 'HD2': '2HD', 'H': 'H', 'HG1': '1HG', 'HG2': '2HG', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB'}, 
             'LEU': {'HD1': 'HD1', 'HD2': 'HD2', 'HD3': 'HD3', 'HD11': '1HD1', 'HD12': '2HD1', 'HD13': '3HD1', 'HD21': '1HD2', 'HD22': '2HD2', 'HD23': '3HD2', 'H': 'H', 'HG': 'HG', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB'}, 
             'ILE': {'HG1': 'HG1', 'HG2': 'HG2','H': 'H', 'HA': 'HA', 'HB': 'HB', 'HG11': '1HG1', 'HG12': '2HG1', 'HG13': '3HG1', 'HG21': '1HG2', 'HG22': '2HG2', 'HG23': '3HG2', 'HD11': '1HD1', 'HD12': '2HD1', 'HD13': '3HD1', 'HD': 'HD1'}, 
             'LYS': {'H': 'H', 'HA': 'HA', 'HB2': '2HB', 'HB1': '1HB', 'HG1': '1HG', 'HG2': '2HG','HD1': '1HD', 'HD2': '2HD','HE1': '1HE', 'HE2': '2HE','HZ1': '1HZ', 'HZ2': '2HZ', 'HZ3': '3HZ'}, 
             'DBS': {'H2': '2H', 'H3': '3H', 'H1': '1H', 'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB': '1HB', 'HB2': '1HB', 'HG1': '1HG', 'HG3': '3HG', 'HG2': '2HG', 'HG': '1HG'},
             'DBR': {'H2': '2H', 'H3': '3H', 'H1': '1H', 'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB': '1HB', 'HB2': '1HB', 'HG1': '1HG', 'HG3': '3HG', 'HG2': '2HG', 'HG': '1HG'}}
#'BDS': {'H2': '2H', 'H3': '3H', 'H1': '1H', 'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB': '1HB', 'HB2': '1HB', 'HG1': '1HG',
#        'HG3': '3HG', 'HG2': '2HG', 'HG': '1HG'},
#'BLR': {'H2': '2H', 'H3': '3H', 'H1': '1H', 'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB': '1HB', 'HB2': '1HB', 'HG1': '1HG',
#        'HG3': '3HG', 'HG2': '2HG', 'HG': '1HG'}}
#this is for the strange HB3 and HB2 naming in nmr structures
terr_dict = {'CYS': {'HB2': '1HB', 'HB3': '2HB'},
             'HIS': {'HE2': '1HE', 'HE3': '2HE', 'HB2': '1HB', 'HB3': '2HB', 'HD2': '1HD', 'HD3': '2HD'}, 
             'GLY': {'HA2': '1HA', 'HA3': '2HA'}, 
             'DAL': {'HB3': '2HB', 'HB2': '1HB'}, 
             'SER': {'HB3': '2HB', 'HB2': '1HB'}, 
             'PRO': {'HD2': '1HD', 'HD3': '2HD', 'HG2': '1HG', 'HG3': '2HG', 'HB2': '1HB', 'HB3': '2HB',}, 
             'ALA': {'H2': '2H', 'H3': '3H', 'H1': '1H', 'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB3': '3HB', 'HB2': '2HB'}, 
             'DHA': {'H': 'H', 'HA': 'HA', 'HB2': '1HB', 'HB3': '2HB'}, 
             'DBU': {'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB3': '3HB', 'HB2': '2HB', 'HG13': '3HG', 'HG12': '2HG', 'HG11': '1HG'}, 
             'DALA': {'H2': '2H', 'H3': '3H', 'H1': '1H', 'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB3': '3HB', 'HB2': '2HB'}, 
             'MET': {'HE1': '1HE', 'HE2': '2HE', 'HE3': '3HE', 'HG2': '1HG', 'HG3': '2HG', 'HB2': '1HB', 'HB3': '2HB'}, 
             'GLU': {'HG2': '1HG', 'HG3': '2HG', 'HB2': '1HB', 'HB3': '2HB'}, 
             'GLN': {'HG2': '1HG', 'HG3': '2HG', 'HB2': '1HB', 'HB3': '2HB', 'HE22': '2HE2', 'HE21': '1HE2'}, 
             'TRP': {'HH2': 'HH2', 'HE1': 'HE1', 'HD1': 'HD1', 'HE3': 'HE3', 'HZ3': 'HZ3', 'HZ2': 'HZ2', 'HA': 'HA', 'HB3': '2HB', 'HB2': '1HB'}, 
             'ASN': {'HD22': '2HD2', 'HD21': '1HD2', 'HB3': '2HB', 'HB1': '1HB'}, 
             'ASP': { 'HB3': '2HB', 'HB2': '1HB'}, 
             'TYR': {'HE2': 'HE2', 'HD2': 'HD2', 'HD1': 'HD1', 'OH': 'OH', 'H': 'H', 'HE1': 'HE1', 'HH': 'HH', 'HA': 'HA', 'HB3': '2HB', 'HB2': '1HB'}, 
             'PHE': {'HE2': 'HE2', 'HD2': 'HD2', 'HD1': 'HD1', 'HZ': 'HZ', 'H': 'H', 'HE1': 'HE1', 'HH': 'HH', 'HA': 'HA', 'HB3': '2HB', 'HB2': '1HB'}, 
             'VAL': {'HG1': 'HG1', 'HG2': 'HG2', 'H': 'H', 'HA': 'HA', 'HB': 'HB', 'HG11': '1HG1', 'HG12': '2HG1', 'HG13': '3HG1', 'HG21': '1HG2', 'HG22': '2HG2', 'HG23': '3HG2'}, 
             'ARG': {'HH21': '1HH2', 'HH22': '2HH2', 'HH11': '1HH1', 'HE12': '2HH1', 'HE': 'HE', 'HD2': '1HD', 'HD3': '2HD', 'HG2': '1HG', 'HG3': '2HG', 'HB3': '2HB', 'HB2': '1HB'}, 
             'LEU': {'HD1': 'HD1', 'HD2': 'HD2', 'HD3': 'HD3', 'HD11': '1HD1', 'HD12': '2HD1', 'HD13': '3HD1', 'HD21': '1HD2', 'HD22': '2HD2', 'HD23': '3HD2', 'HG': 'HG', 'HB3': '2HB', 'HB2': '1HB'}, 
             'ILE': {'HB': 'HB', 'HG12': '1HG1', 'HG13': '2HG1', 'HG21': '1HG2', 'HG22': '2HG2', 'HG23': '3HG2', 'HD11': '1HD1', 'HD12': '2HD1', 'HD13': '3HD1', 'HD': 'HD1'}, 
             'LYS': {'HB3': '2HB', 'HB2': '1HB', 'HG2': '1HG', 'HG3': '2HG','HD2': '1HD', 'HD3': '2HD','HE2': '1HE', 'HE3': '2HE','HZ1': '1HZ', 'HZ2': '2HZ', 'HZ3': '3HZ'}, 
             'DBS': {'H2': '2H', 'H3': '3H', 'H1': '1H', 'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB': '1HB', 'HB2': '1HB', 'HG1': '1HG', 'HG3': '3HG', 'HG2': '2HG', 'HG': '1HG'},
             'DBR': {'H2': '2H', 'H3': '3H', 'H1': '1H', 'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB': '1HB', 'HB2': '1HB', 'HG1': '1HG', 'HG3': '3HG', 'HG2': '2HG', 'HG': '1HG', 'HG21': '1HG', 'HG22': '2HG', 'HG23': '3HG'}}
#'BDS': {'H2': '2H', 'H3': '3H', 'H1': '1H', 'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB': '1HB', 'HB2': '1HB', 'HG1': '1HG',
#        'HG3': '3HG', 'HG2': '2HG', 'HG': '1HG'},
#'BLR': {'H2': '2H', 'H3': '3H', 'H1': '1H', 'H': 'H', 'HA': 'HA', 'HB1': '1HB', 'HB': '1HB', 'HB2': '1HB', 'HG1': '1HG',
#        'HG3': '3HG', 'HG2': '2HG', 'HG': '1HG'}}
#temporarily saying HG = 1HG for BDS, also has 1HB, not HB or 2HB
#H HA HH

out_str = ''
filetype = 'norm'
if argument.name != 'norm':
    filetype = 'terrible'
    print('using odd naming convention')
count = 0
loose = bool(argument.loose)
start = 0

seq_dict = {}
with open(argument.seq, 'r') as infile:
    first_line = infile.readline()
    split = first_line.split()
    count = 0
    if split[0] == '2KT':
        start = -1
    for resi in split:
        count += 1
        seq_dict[str(count + start)] = resi

def make_constraint(res1, res2, atom1, atom2, center, lower, upper):
    #print(res1, res2, atom1, atom2, seq_dict[res1], seq_dict[res2])
    if atom1 == 'HBR':
        atom1 = 'HB3'
    if atom2 == 'HBR':
        atom2 = 'HB3'
    if atom1 == 'HBS':
        atom1 = 'HB2'
    if atom2 == 'HBS':
        atom2 = 'HB2'
    if '*' in atom1:
        if 'DB' not in seq_dict[res1]:
            atom1 = 'Q' + atom1.replace('*', '')[1:]
        else:
            atom1 = atom1.replace('*', '')
    if '*' in atom2:
        if 'DB' not in seq_dict[res2]:
            atom2 = 'Q' + atom2.replace('*', '')[1:]
        else:
            atom2 = atom2.replace('*', '')
    ambig = False
    if 'Q' in atom1 or 'Q' in atom2: #if ambig, Rosetta should be fine
        ambig = True
    if 'Q' not in atom1 and '#' not in atom1 and atom1[-1] in ['1', '2', '3']:
        #print(seq_dict[res1], atom1, res1)
        if filetype == 'norm':
            #atom1 = atom1[-1] + atom1[:-1]
            atom1 = name_dict[seq_dict[res1]][atom1]
        else:
            atom1 = terr_dict[seq_dict[res1]][atom1]
    if 'Q' not in atom2 and '#' not in atom2 and atom2[-1] in ['1', '2', '3']:
        #print(seq_dict[res2], res2, atom2)
        if filetype == 'norm':
            #atom2 = atom2[-1] + atom2[:-1]
            atom2 = name_dict[seq_dict[res2]][atom2]
        else:
            if atom2 not in terr_dict[seq_dict[res2]].keys():
                print('missing key for: ')
                print(seq_dict[res2], res2, atom2)
            atom2 = terr_dict[seq_dict[res2]][atom2]
    if '#' in atom1:
        atom1 = atom1.replace('#', '') #rosetta can handle HB or QB main/source/src/core/scoring/constraints/AmbiguousNMRDistanceConstraint.cc
        ambig = True
    if '#' in atom2:
        atom2 = atom2.replace('#', '')
        ambig = True
    if seq_dict[res1] in ['DBR', 'DBS'] and atom1 == 'HB':
        atom1 = '1HB'
    if seq_dict[res2] in ['DBR', 'DBS'] and atom2 == 'HB':
        atom2 = '1HB'
    if atom1 == 'HN':
        atom1 = 'H'
    if atom2 == 'HN':
        atom2 = 'H'
    if [atom1, seq_dict[res1]] in [['HD', "ILE"], ['HG', 'DBS'], ['HB', 'DBS']]:
        #== 'HD' and res1 == '13':
        #print(line)
        #print(name_dict[seq_dict[res1]][atom1], seq_dict[res1], atom1, res1)
        atom1 = name_dict[seq_dict[res1]][atom1]
    if [atom2, seq_dict[res2]] in [['HD', "ILE"], ['HG', 'DBS'], ['HB', 'DBS']]:
        atom2 = name_dict[seq_dict[res2]][atom2]
    if loose:
        #increasing max distance by 1 angsrom and lowering standard deviation value per Jens
        if ambig:
            #count += 1
            return 'AmbiguousNMRDistance '+atom1+' '+res1+' '+atom2+' '+res2+' BOUNDED '+str(0.0)+' '+str(round(center + upper + 1.0, 2))+' 0.5 0.5 count\n'
        else:
            #count += 1
            return 'AtomPair '+atom1+' '+res1+' '+atom2+' '+res2+' BOUNDED '+str(0.0)+' '+str(round(center + upper + 1.0, 2))+' 0.5 0.5 count\n'
    else:
        if ambig:
            #count += 1
            return 'AmbiguousNMRDistance '+atom1+' '+res1+' '+atom2+' '+res2+' BOUNDED '+str(center-lower)+' '+str(round(center + upper, 2))+' 1.0 0.5 count\n'
        else:
            #count += 1
            return 'AtomPair '+atom1+' '+res1+' '+atom2+' '+res2+' BOUNDED '+str(center-lower)+' '+str(round(center + upper, 2))+' 1.0 0.5 count\n'

def fix_naming(line):
    split = line.split()
    line = line.replace('HB3', 'HB2').replace('HB2', 'HB1').replace('HB1', 'HB2')
    line = line.replace('HG3', 'HG2').replace('HG2', 'HG1').replace('HG1', 'HG2')
    line = line.replace('HD3', 'HD2').replace('HD2', 'HD1').replace('HD1', 'HD2')
    if 'DBB' in line:
        if any(len(x) == 4 and 'HG2' in x for x in split):
            line = line.replace('HG2', 'HG')
    return line
    

aa_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
           'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 
           'DBB', 'DBU', 'DHA', 'DBS', 'DBR', 'DAL', 'ABA'] #'BDS', 'BLR'

with open(argument.input, 'r') as infile:
    for line in infile:
        if 'HA3' in line or 'HD3' in line:
            filetype = 'terrible' #some files use HB3 and HB2 instead of HB1 and HB2
            print('using odd naming convention')
            break

conn_dict = {}
if argument.conn != '':
    with open(argument.conn, 'r') as infile:
        for raw_line in infile:
            line = raw_line.strip()
            split = line.split(',')
            #non-cys residue always first
            conn_dict[int(split[0])] = split[1]
print(conn_dict)

def fix_conn_naming(res1, res2, atom1, atom2):
    if int(res1) in conn_dict.keys():
        # print(res1)
        if conn_dict[int(res1)][0:2] == 'DB':  # methyllanthionine, based on Rosetta naming and patch
            if 'HB' in atom1:
                atom1 = 'HB1'
        else:  # ALA or DALA, NMR NOE assignment name and structure name differ (2KTO)
            if atom1 == 'HB2':
                atom1 = 'HB1'
            elif atom1 == 'HB3':
                atom1 = 'HB2'
    if int(res2) in conn_dict.keys():
        # print(res2)
        if conn_dict[int(res2)][0:2] == 'DB':  # methyllanthionine
            if 'HB' in atom2:
                atom2 = 'HB1'
        else:  # ALA or DALA, NMR NOE assignment name and structure name differ (2KTO)
            if atom2 == 'HB2':
                atom2 = 'HB1'
            elif atom2 == 'HB3':
                atom2 = 'HB2'
    return atom1, atom2


with open (argument.input, 'r') as infile:
    count = 0
    for raw_line in infile:
        line = raw_line.strip()
        if argument.fix:
            line = fix_naming(line)
        split = line.split()
        if len(line) == 0 or line[0] == '#':
            pass
        elif line[0:6] == 'assign': #this is for XPLOR type files
            res1 = str(int(split[2]) + start)
            res2 = str(int(split[7]) + start)
            atom1 = split[5][:-1]
            atom2 = split[10][:-1]
            center = float(split[11])
            lower = float(split[12])
            upper = float(split[13])
            #ambig = False
            #print(line)
            atom1, atom2 = fix_conn_naming(res1, res2, atom1, atom2)
            if center > 2.3 :
                temp_str = make_constraint(res1, res2, atom1, atom2, center, lower, upper)
                count += 1
                out_str += temp_str.replace('count', str(count))
        elif split[1] in aa_list and len(split) <= 13: #this is for DYANA type files
            res1 = str(int(split[0]) + start)
            res2 = str(int(split[3]) + start)
            atom1 = split[2]
            atom2 = split[5]
            center = float(split[6])
            lower = center - 1.8
            upper = 0.0
        #    temp_str = make_constraint(res1, res2, atom1, atom2, center, lower, upper)
        #    count += 1
        #    out_str += temp_str.replace('count', str(count))
            atom1, atom2 = fix_conn_naming(res1, res2, atom1, atom2)
            if center > 2.3 :
                temp_str = make_constraint(res1, res2, atom1, atom2, center, lower, upper)
                count += 1
                out_str += temp_str.replace('count', str(count))
        elif len(split) == 64 and split[8] in aa_list and split[18] in aa_list: #2KTO
            res1 = str(int(split[7]) + start)
            res2 = str(int(split[17]) + start)
            atom1 = split[51]
            atom2 = split[59]
            center = float(split[28])
            lower = center - 1.8
            upper = 0.0
            #lets fix any connection naming
            atom1, atom2 = fix_conn_naming(res1, res2, atom1, atom2)
            #if int(res1) in conn_dict.keys():
            #    #print(res1)
            #    if conn_dict[int(res1)][0:2] == 'DB': #methyllanthionine, based on Rosetta naming and patch
            #        if 'HB' in atom1:
            #            atom1 = 'HB1'
            #    else: #ALA or DALA, NMR NOE assignment name and structure name differ (2KTO)
            #        if atom1 == 'HB2':
            #            atom1 = 'HB1'
            #        elif atom1 == 'HB3':
            #            atom1 = 'HB2'
            #if int(res2) in conn_dict.keys():
            #    #print(res2)
            #    if conn_dict[int(res2)][0:2] == 'DB': #methyllanthionine
            #        if 'HB' in atom2:
            #            atom2 = 'HB1'
            #    else: #ALA or DALA, NMR NOE assignment name and structure name differ (2KTO)
            #        if atom2 == 'HB2':
            #            atom2 = 'HB1'
            #        elif atom2 == 'HB3':
            #            atom2 = 'HB2'
            if center > 2.3 :
                temp_str = make_constraint(res1, res2, atom1, atom2, center, lower, upper)
                count += 1
                out_str += temp_str.replace('count', str(count))
        else:
            #for i in range(len(split)):
            #    print(i, split[i])
            #break
            print(f'{split[8]} or {split[18]} not an AA')

basename = argument.input.split('.')[0]

with open(basename+'.cst', 'w') as outfile:
    outfile.write(out_str[:-1])
