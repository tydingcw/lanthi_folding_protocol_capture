from rdkit import Chem
import argparse
#import numpy as np

parser = argparse.ArgumentParser(description = 'take in a pdb file and rename lanthionine residues based on steriochemistry')
parser.add_argument('-i', '--input', help = 'input pdb', required = True)
parser.add_argument('-c', '--chain', help = 'chain of interest', default='A')
parser.add_argument('-p', '--prefix', help = 'prefix for output', default='')
parser.add_argument('-w', '--write', help = 'write connect file?', default=True)
parser.add_argument('-k', '--keep', help = 'keep other chains in output?', default=False)
parser.add_argument('-f', '--file', help = 'pdb file to define ring topology instead of using input', default='')
#parser.add_argument('-a', '--alphas', help = 'prepend ca carbons', default=0, type=int)
argument=parser.parse_args()

# def get_dihedral(atom1, atom2, atom3, atom4):
#     #https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
#     """Praxeolitic formula
#     1 sqrt, 1 cross product"""
#     p0 = np.array(atom1)
#     p1 = np.array(atom2)
#     p2 = np.array(atom3)
#     p3 = np.array(atom4)

#     b0 = -1.0*(p1 - p0)
#     b1 = p2 - p1
#     b2 = p3 - p2

#     # normalize b1 so that it does not influence magnitude of vector
#     # rejections that come next
#     b1 /= np.linalg.norm(b1)

#     # vector rejections
#     # v = projection of b0 onto plane perpendicular to b1
#     #   = b0 minus component that aligns with b1
#     # w = projection of b2 onto plane perpendicular to b1
#     #   = b2 minus component that aligns with b1
#     v = b0 - np.dot(b0, b1)*b1
#     w = b2 - np.dot(b2, b1)*b1

#     # angle between v and w in a plane is the torsion angle
#     # v and w may not be normalized but that's fine since tan is y/x
#     x = np.dot(v, w)
#     y = np.dot(np.cross(b1, v), w)
#     angle = np.degrees(np.arctan2(y, x))
#     return angle

ripp_file = argument.input
ripp_chain = argument.chain
if argument.file != '':
    ripp_file =  argument.file
    ripp_chain = 'A'
    

ripp = Chem.MolFromPDBFile(ripp_file)
#ripp = Chem.RWMol(ripp)
#ripp.AddBond(31, 70)
#ripp.AddBond(111, 164)
#ripp.AddBond(19, 42)
#ripp.AddBond(61, 90)
#ripp = Chem.rdchem.Mol(ripp)
#print(type(ripp))
ca_list = []
for i in range(len(ripp.GetAtoms())):
    if 'CA' in ripp.GetAtomWithIdx(i).GetPDBResidueInfo().GetName() and ripp.GetAtomWithIdx(i).GetPDBResidueInfo().GetChainId() == ripp_chain:
        #print(ripp.GetAtomWithIdx(i).GetPDBResidueInfo().GetChainId())
        #print(ripp.GetAtomWithIdx(i).GetPDBResidueInfo().GetName())
        #print(i)
        #print(ripp.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName())
        ca_list.append(i)
#print(ca_list)
cb_list = []
for i in range(len(ripp.GetAtoms())):
    if 'CB' in ripp.GetAtomWithIdx(i).GetPDBResidueInfo().GetName() and ripp.GetAtomWithIdx(i).GetPDBResidueInfo().GetChainId() == ripp_chain:
        #print(ripp.GetAtomWithIdx(i).GetPDBResidueInfo().GetName())
        #print(i)
        cb_list.append(i)
#print(cb_list)
#methio = ripp.GetSubstructMatches(Chem.MolFromSmiles('CCSC(C)C'))
#print(ripp.GetSubstructMatches(Chem.MolFromSmiles('CS')))
methio = ripp.GetSubstructMatches(Chem.MolFromSmiles('NCCSCCN'))
#print(methio)
#cystine should be R
methio_ca = []
methio_ca_dict = {}
for item in methio:
    for index in item:
        if index in ca_list:
            methio_ca.append(index)
            if len(methio_ca) % 2 == 0:
                methio_ca_dict[methio_ca[-1]] = methio_ca[-2]
                methio_ca_dict[methio_ca[-2]] = methio_ca[-1]
#print(methio_ca)
#print(methio_ca_dict)
methio = ripp.GetSubstructMatches(Chem.MolFromSmiles('SC(C)CN'))
#print(methio)
methio_cb = []
#methio_cb_dict = {}
for item in methio:
    for index in item:
        if index in cb_list:
            methio_cb.append(index)
            #if len(methio_cb) % 2 == 0:
            #    methio_cb_dict[methio_cb[-1]] = methio_cb[-2]
            #    methio_cb_dict[methio_cb[-2]] = methio_cb[-1]
#print(methio_cb)
#print(Chem.FindMolChiralCenters(ripp, includeUnassigned=True, useLegacyImplementation=False))
#print(ca_list)
#print(cb_list)

#first steriocenter is for D/L second is for R/S on methylanthionine
#R=L becuase of the sulfur in the lanthionine ring
og_name_dict={('R', 'S'): 'BLS', ('S', 'R'): 'BDR', ('R', 'R'): 'BLR', ('S', 'S'): 'BDS'}
#will let rosetta handle the D/L
name_dict={('R', 'S'): 'DBS', ('S', 'R'): 'DBR', ('R', 'R'): 'DBR', ('S', 'S'): 'DBS'}
connect_dict = {}
pdb_dict = {}

chiral_list = Chem.FindMolChiralCenters(ripp, includeUnassigned=True, useLegacyImplementation=False)
for item in chiral_list:
    ind = item[0]
    if ind in ca_list:
        if ind in methio_ca:
            pdb_ind = ca_list.index(ind)+1
            pdb_dict[ind]=pdb_ind

lanthionine_list = []
methyllanthionine_list = []
for item in chiral_list:
    ind = item[0]
    if ind in ca_list:
        if ind in methio_ca:
            pdb_ind = ca_list.index(ind) + 1 #+ argument.alphas
            #print(pdb_ind, item, 'methionine ca')
            #print(ripp.GetAtomWithIdx(ind).GetPDBResidueInfo().GetResidueName())
            connect_dict[pdb_ind] = {}
            connect_dict[pdb_ind]['CA'] = item
            connect_dict[pdb_ind]['connect']=pdb_dict[methio_ca_dict[ind]]
            connect_dict[pdb_ind]['CB'] = None
            connect_dict[pdb_ind]['name'] = ripp.GetAtomWithIdx(ind).GetPDBResidueInfo().GetResidueName()
            
            if connect_dict[pdb_ind]['name'] in ['DBB', 'DAL', 'ALA']: #or connect_dict[pdb_ind]['name'] == 'DAL':
                lanthionine_list.append(pdb_ind)
                if connect_dict[pdb_ind]['name'] == 'DBB':
                    methyllanthionine_list.append(pdb_ind)
            
        else:
            pass
            #print(ca_list.index(ind)+1, item)
            #print(ripp.GetAtomWithIdx(ind).GetPDBResidueInfo().GetResidueName())
    elif ind in cb_list:
        if ind in methio_cb:
            pdb_ind = ca_list.index(ind-3)+1
            #print(ca_list.index(ind-3)+1, item, 'methionine cb')
            #print(ripp.GetAtomWithIdx(ind).GetPDBResidueInfo().GetResidueName())
            connect_dict[pdb_ind]['CB'] = item

#print(connect_dict)
#print(lanthionine_list)
#print(methio_ca)
#print(methio_ca_dict)

#determine E/Z for DBU
#take dihedral with coords from N, CA, CB, CG
#dbu_dict = {}
#dbu_ind = 0

# with open(argument.input, 'r') as in_file:
#     for raw_line in in_file:
#         line=raw_line.strip()
#         split = line.split()
#         if line[0:6] == 'HETATM' and split[4] == argument.chain and split[3] == 'DBU':
#             if split[5] not in dbu_dict.keys():
#                 dbu_dict[split[5]] = {}
#             if split[2] == 'N':
#                 dbu_dict[split[5]]['N'] = [float(split[6]), float(split[7]), float(split[8])]
#             elif split[2] == 'CA':
#                 dbu_dict[split[5]]['CA'] = [float(split[6]), float(split[7]), float(split[8])]
#             elif split[2] == 'CB':
#                 dbu_dict[split[5]]['CB'] = [float(split[6]), float(split[7]), float(split[8])]
#             elif split[2] == 'CG':
#                 dbu_dict[split[5]]['CG'] = [float(split[6]), float(split[7]), float(split[8])]
#             tmp_keys = dbu_dict[split[5]].keys()
#             if 'conf' not in tmp_keys and 'N' in tmp_keys and 'CA' in tmp_keys and 'CB' in tmp_keys and 'CG' in tmp_keys:
#                 #get_dihedral(atom1, atom2, atom3, atom4)
#                 sub_dict = dbu_dict[split[5]]
#                 dih = get_dihedral(sub_dict['N'], sub_dict['CA'], sub_dict['CB'], sub_dict['CG'])
#                 #print(dih)
#                 if -90 < dih < 90:
#                     dbu_dict[split[5]]['conf'] = 'Z'
#                 else:
#                     dbu_dict[split[5]]['conf'] = 'E'

#output a file specifying the connections for lanthionine/methylanthionine
#output a file with DBB renamed based on steriochemistry
connect_str = ''
chiral_str = ''
for ind in lanthionine_list:
    resname=connect_dict[ind]['name']
    connect_ind=connect_dict[ind]['connect']
    connect_name=connect_dict[connect_ind]['name']
    if resname == 'DBB':
        resname = name_dict[(connect_dict[ind]['CA'][1], connect_dict[ind]['CB'][1])]
    #print(f'{ind},{name},CONNECT,{connect_ind},{connect_name}')
    connect_str += f'{ind},{resname},CONNECT,{connect_ind},{connect_name}\n'
    chiral_str += f'{ind},{resname},CONNECT,{connect_ind},{connect_name}\n'

out_str = ''
link_str = ''
nmr=False
model=1
split_input = argument.input.split('.')

#example of link in pdb
#LINK         CB  DBB A  10                 SG  CYS A  32     1555   1555  1.82
#HETATM  101  HB2 BDS A   7      10.739  -0.423  -2.909  1.00  0.00           H

supported = True

ca_init = -1 #using this to catch when residue numbering does not match ca_index
with open(argument.input, 'r') as in_file:
    for raw_line in in_file:
        #print(ca_init, 'ca')
        line=raw_line.strip()
        if ca_init == -1 and (line[0:4] == "ATOM" or line[0:6] == 'HETATM') and 'N' in line: # and 'CA' in line:
            ca_init += int(line.split()[5])
            #print('add to ca_init')
            #print(line)
            #print(line.split()[5])
        if line[0:4] == 'LINK' and 'CB' in line and 'SG' in line:
            link_str += raw_line
            link_str = link_str.replace("DBB", resname)
        elif line[0:5] == 'MODEL':
            nmr=True
            model = line.split()[1]
            out_str = ''
        elif (line[0:4] == "ATOM" or line[0:6] == 'HETATM') and (line.split()[4] == argument.chain or argument.keep):
            if line.split()[3] == 'DBB': #['17:20']
                resn = int(line.split()[5]) - ca_init
                #print('breakpoint')
                #print(ca_init)
                #print(name_dict)
                #print(line)
                #print(connect_dict)
                #print(resn)
                resname = name_dict[(connect_dict[resn]['CA'][1], connect_dict[resn]['CB'][1])]
                add_line = line[:17] + resname + line[20:] + '\n'
                if 'HB' in add_line:
                    ind = add_line.index('HB')
                    add_line = add_line[:ind] + 'HB1' + add_line[ind+3:]
                out_str += add_line
            #elif line.split()[3] == 'DBU':
            #    #print(dbu_dict)
            #    if dbu_dict[line.split()[5]]['conf'] == 'Z':
            #        out_str += line.replace('DBU', 'ZBU') + '\n'
            #    elif dbu_dict[line.split()[5]]['conf'] == 'E':
            #        out_str += line.replace('DBU', 'EBU') + '\n'
            #    else:
            #        print('error dih is: '+dbu_dict[line.split()[5]]['conf'])
            elif line[0:4] == "ATOM":
                out_str += line + '\n'
            elif line[0:6] == 'HETATM' and line.split()[3] in ['DHA', 'DAL', 'DBU']:
                out_str += line + '\n'
            elif line[0:6] == 'HETATM' and line.split()[3] in ['2KT']:
                pass #delete for now
            else: 
                #not an CAA, DBB, DAL, DHA, DBB
                print('NCAA that is not one of DBB, DAL, DHA, DBU was detected ('+ line.split()[3] +'), so breaking')
                supported = False
                break
        elif line[0:3] == 'TER' and out_str != '':
            if nmr:
                out_name = f'{split_input[0]}_{model}_rename.pdb'
                with open(argument.prefix + out_name, 'w') as out_file:
                    out_file.write(link_str+out_str)
            else:
                out_name = f'{split_input[0]}_rename.pdb'
                with open(argument.prefix + out_name, 'w') as out_file:
                    out_file.write(link_str+out_str)
                    
if supported and type(argument.write)==bool:
    with open(argument.prefix + split_input[0]+'_conn.txt', 'w') as out_file:
        out_file.write(connect_str[:-1])
    with open(argument.prefix + split_input[0]+'_chiral.txt', 'w') as out_file:
        out_file.write(chiral_str[:-1])
