#this script is used to generate gaussian input and parse results
from argparse import ArgumentParser
import numpy as np

class Rotamer:
    def __init__(self):
        self.atom_list = []
        self.bond_dict = {}
        #self.sulfur_pos = -1
        self.conf_score = None
        self.clash_score = None
    
    def add_atom(self, atom):
        self.atom_list.append(atom)
        #if atom.name == 'S':
        #    self.sulfur_pos = len(self.atom_list)-1
    
    def add_bond(self, atom1, atom2):
        atom1 = int(atom1)
        atom2 = int(atom2)
        if atom1-1 in self.bond_dict.keys():
            self.bond_dict[atom1-1].append(atom2-1)
        else:
            self.bond_dict[atom1-1] = [atom2-1]
        if atom2-1 in self.bond_dict.keys():
            self.bond_dict[atom2-1].append(atom1-1)
        else:
            self.bond_dict[atom2-1] = [atom1-1]

class Atom:
    def __init__(self, name, x, y, z):
        self.name = name
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
    
    def get_xyz(self):
        return [self.x, self.y, self.z]
    
def main():
    help_text = "Must specify one of: gjf - generate gjf from sdf; "
    help_text += "freq - check if freqs are all positive; "
    help_text += "extract - extract gout coords from log file; "
    help_text += "zmat - convert to zmat gjf (required -i *.sdf -d 1,2,3,4) "
    #add functionality to generate partial charge file
    #add gen_zmat code for dih,angle,bond scans
    #add options to specify theory and if freq is desired
    parser = ArgumentParser(description="command line arguments")
    parser.add_argument('-i', '--input', required=True, help="Input file")
    parser.add_argument('-m', '--method', required=True, help=help_text)
    parser.add_argument('-o', '--output', default='default', help="Path for output if needed")
    parser.add_argument('-p', '--processors', default='12', help="Number of processors, default is 12")
    parser.add_argument('-r', '--ram', default='8000MB', help="Amount of ram, default is 8000MB")
    #parser.add_argument('-c', '--connect', default='def', help="sdf file to define connectivity of xyz")
    parser.add_argument('-d', '--dihedral', default='def', help="atom ids that define the dihedral of interest ex: 1,2,3,4")
    args=parser.parse_args()
    
    if args.method == 'gjf':
        if args.input[-4:] == '.sdf':
            rot_list = make_rotamer_objs(args.input)
            rotamer_to_gjf(rot_list[0], args.output.split('.')[0], args.ram, args.processors)
        else:
            print('for method gjf, a .sdf file was not provided')
    elif args.method == 'freq':
        if args.input[-4:] == '.log':
            freq_check(args.input)
        else:
            print('for method freq, a .log file was not provided')
    elif args.method == 'extract':
        if args.input[-4:] == '.log':
            extract_gout_coordinate(args.input, args.output.split('.')[0])
        else:
            print('for method extract, a .log file was not provided')
    elif args.method =='zmat':
        #input file should be xyz of optimized geometry
        #need to input an sdf as well to get the connectivity for zmat
        #need to specify dihedral of interest as well or whatnot
        #if args.connect == 'def':
        #    print('when constructing zmatrix, need to specify sdf file for connectivity (flag -c)')
        if args.input[-4:] != '.sdf':
            print('sdf not provided')
        if args.dihedral == 'def':
            print('when constructing zmatrix, need to specify the initial dihedral (flag -d)')
        elif args.dihedral.count(',') != 3 or args.dihedral.count(' ') != 0:
            print('-d flag should have four numbers separated by three , and no spaces')
        if args.input[-4:] == '.sdf' and args.dihedral != 'def' and args.dihedral.count(',') == 3 or args.dihedral.count(' ') == 0:
            #rot_connect = make_rotamer_objs(args.connect)
            #rot_connect = rot_connect[0]
            rot_zmat = make_rotamer_objs(args.input)[0]
            #print(args.input)
            #print(rot_zmat)
            #rot_zmat = rot_zmat[0]
            #flag = False
            #if len(rot_connect.atom_list) == len(rot_zmat.atom_list):
            #    for i in range(len(rot_connect.atom_list)):
            #        if rot_connect.atom_list[i].name != rot_zmat.atom_list[i].name:
            #            flag = True
            #            print('at pos '+str(i)+' connect has '+rot_connect.atom_list[i].name+'and zmat has '+rot_zmat.atom_list[i].name)
            #else:
            #    print('unequal atom list len')
            #    flag = True
            #if flag:
            #    print('atom lists do not match')
            #else:
            #rot_zmat.bond_dict = rot_connect.bond_dict
            split = args.dihedral.split(',')
            dihedral=[int(x) for x in split]
            gen_zmat(rot_zmat, dihedral, args.output.split('.')[0], set_dih = 0.0, set_steps = 35, set_change = 10.0, 
                         nproc=args.processors, theory='M052X/6-311+G(d,p)', orig_dih=True)
            print('ran gen_zmat')
    else:
        print("Bad input for method (-m). " + help_text)

def make_rotamer_objs(filename):
    rotamer_list = []
    read_frame = False
    get_conf_score = False
    get_clash_score = False
    with open(filename, 'r') as file:
        print(filename)
        for raw_line in file:
            line = raw_line.strip()
            line_list = line.split()
            if len(line_list)>=11 and not read_frame:
                read_frame = True
                rotamer_list.append(Rotamer())
            elif 'END' in raw_line:
                read_frame = False
            elif read_frame and len(line_list) == 16:
                #have an atom
                new_atom = Atom(line_list[3], line_list[0], 
                               line_list[1], line_list[2])
                rotamer_list[len(rotamer_list)-1].add_atom(new_atom)
            elif read_frame and len(line_list) == 7:
                #have a bond
                rotamer_list[len(rotamer_list)-1].add_bond(line_list[0], 
                                                           line_list[1])
            if get_conf_score:
                get_conf_score = False
                rotamer_list[len(rotamer_list)-1].conf_score = float(line)
            elif '<ConfScore>' in line:
                get_conf_score = True
            
            if get_clash_score:
                get_clash_score = False
                rotamer_list[len(rotamer_list)-1].clash_score = float(line)
            elif '<ConfClashScore>' in line:
                get_clash_score = True
            
        return rotamer_list

#make gjf
#https://comp.chem.umn.edu/info/dft.htm
def rotamer_to_gjf(rot, path, ram, proc):
    out_text = ''
    new_line = '\n'
    #freq
    route = '%chk='+path+'.chk'+new_line+'%mem='+ram+new_line+'%nprocshared='+proc+new_line+'# opt M052X/6-311+G(d,p) scrf=smd'
    title = 'Title Card Required'
    chrgspin = '0 1'
    out_text = route + new_line + new_line + title + new_line + new_line + chrgspin + new_line
    for atom in rot.atom_list:
        line = atom.name + ' ' + str(atom.x) + ' ' + str(atom.y) + ' ' + str(atom.z) + new_line
        out_text = out_text + line
    out_text = out_text + new_line
    
    with open(path+'.gjf', 'w') as file:
        file.write(out_text)
    
    #return out_text

def freq_check(gout_path):
    all_freq = []
    with open(gout_path, 'r') as file:
        for raw_line in file:
            line = raw_line.strip()
            if 'Frequencies' in line:
                split=line.split()
                for x in split[2:]:
                    all_freq.append(x) 
    if all(all_freq) > 0:
        print('All Freq Positive :)')
    else:
        print('Negative Freq Found!!!')

def extract_gout_coordinate(gout_path, filename=''):
    out_rot=Rotamer()
    #coords = ''
    out_string=''
    with open(gout_path, 'r') as raw_file:
        file = raw_file.readlines()
        read_open = -1
        for i in range(len(file)):
            if 'Standard orientation:' in file[i]:
                read_open = i
        read_close = -1
        for j in range (read_open + 5, len(file)):
            if '-----' in file[j]:
                read_close = j
                break
        atom_num_dict = {'1': 'H', '6': 'C', '7': 'N', '8': 'O', '9': 'F', '16': 'S', '17': 'Cl'}
        #out_text = ''
        for line in file[read_open+5 : read_close] :
            split = line.strip().split()
            #print(atom_num_dict[split[1]]+' '+split[3]+' '+split[4]+' '+split[5])
            #coords = coords + atom_num_dict[split[1]]+' '+split[3]+' '+split[4]+' '+split[5]+'\n'
            out_rot.add_atom(Atom(atom_num_dict[split[1]], split[3], split[4], split[5]))
            if filename != '':
                out_string += atom_num_dict[split[1]]+' '+split[3]+' '+split[4]+' '+split[5]+'\n'
    #pass
    print('Gaussian coordinates extracted')
    if filename != '':
        #add check to see if xyz or not
        with open(filename + '.xyz', 'w') as out_file:
            out_file.write(out_string)
    #return coords
    return out_rot

def get_dist(atom1, atom2):
    dist = ((atom1[0]-atom2[0])**2 + (atom1[1]-atom2[1])**2 + 
            (atom1[2]-atom2[2])**2)**(0.5)
    return dist

def get_angle(atom1, atom2, atom3):
    #https://stackoverflow.com/questions/35176451/python-code-to-calculate-angle-between-three-point-using-their-3d-coordinates
    a = np.array(atom1)
    b = np.array(atom2)
    c = np.array(atom3)
    
    #print(a)
    #print(b)
    #print(c)

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.degrees(np.arccos(cosine_angle))
    
    return angle

def get_dihedral(atom1, atom2, atom3, atom4):
    #https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = np.array(atom1)
    p1 = np.array(atom2)
    p2 = np.array(atom3)
    p3 = np.array(atom4)

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    angle = np.degrees(np.arctan2(y, x))
    return angle

def gen_zmat(rotamer, dihedral, path, set_dih = 0.0, set_steps = 36, set_change = 10.0, 
             nproc='12', theory='#M062X/def2TZVP', orig_dih=True):
    #use get_dist, get_angle, get_dihedral
    defined_atoms=[]
    tot_atoms = len(rotamer.atom_list)
    undefined_atoms = [x for x in range(tot_atoms)]
    new_line = '\n'
    out_string = ''
    z_name={}
    z_atom_tot={}
    for atom in rotamer.atom_list: z_atom_tot[atom.name]=0;
    
    #intitialize names
    for ind in undefined_atoms:
        atom_type = rotamer.atom_list[ind].name
        z_atom_tot[atom_type] += 1
        z_name[ind]=atom_type+str(z_atom_tot[atom_type])
    
    #find four connected atoms
    '''start_chains = []
    def init_chain(start, chain_len=4, past_list = [0]):
        past = [x for x in past_list]
        if len(rotamer.bond_dict[start]) > 1 and len(past) < chain_len:
            for i in rotamer.bond_dict[start]:
                #find the next chain element, make sure not to repeat
                if i not in past:
                    temp_past = copy.copy(past)
                    temp_past.append(i)
                    #print(str(i)+' '+str(temp_past))
                    if len(temp_past) != 4:
                        init_chain(i, chain_len, temp_past)
                    else:
                        #these are the dihedrals of interst
                        #print(temp_past)
                        start_chains.append(temp_past)
                        #print(start_chains)
                else:
                    #we have visited here already, don't add to chain
                    pass
        if len(rotamer.bond_dict[start]) == 1:
            #this is a terminus, so stop chain
            pass
            return -1'''
    
    dih_start = [x-1 for x in dihedral]
    #first atom
    out_string = z_name[dih_start[0]]+new_line
    defined_atoms.append(dih_start[0])
    #second atom
    out_string += z_name[dih_start[1]]+' '+z_name[dih_start[0]]+' '+str(get_dist(rotamer.atom_list[dih_start[0]].get_xyz(), 
                                                              rotamer.atom_list[dih_start[1]].get_xyz()))+new_line
    defined_atoms.append(dih_start[1])
    #third atom
    out_string += z_name[dih_start[2]]+' '+z_name[dih_start[1]]+' '+str(get_dist(rotamer.atom_list[dih_start[2]].get_xyz(), 
                                                              rotamer.atom_list[dih_start[1]].get_xyz()))+' '
    out_string += z_name[dih_start[0]]+' '+str(get_angle(rotamer.atom_list[dih_start[2]].get_xyz(), 
                                                         rotamer.atom_list[dih_start[1]].get_xyz(), 
                                                         rotamer.atom_list[dih_start[0]].get_xyz()))+new_line
    defined_atoms.append(dih_start[2])
    #fourth atom
    out_string += z_name[dih_start[3]]+' '+z_name[dih_start[2]]+' '+str(get_dist(rotamer.atom_list[dih_start[3]].get_xyz(), 
                                                              rotamer.atom_list[dih_start[2]].get_xyz()))+' '
    out_string += z_name[dih_start[1]]+' '+str(get_angle(rotamer.atom_list[dih_start[3]].get_xyz(), 
                                                         rotamer.atom_list[dih_start[2]].get_xyz(), 
                                                         rotamer.atom_list[dih_start[1]].get_xyz()))+' '
    #choosing to use the original dihedral in zmat, or reset it
    if orig_dih:
        out_string += z_name[dih_start[0]]+' '+str(get_dihedral(rotamer.atom_list[dih_start[3]].get_xyz(), 
                                                                rotamer.atom_list[dih_start[2]].get_xyz(), 
                                                                rotamer.atom_list[dih_start[1]].get_xyz(), 
                                                                rotamer.atom_list[dih_start[0]].get_xyz()))
        out_string += new_line
    else:
        #if relax:
        out_string += z_name[dih_start[0]]+' '+str(set_dih)+new_line
        #else:
        #    out_string += z_name[dih_start[0]]+' '+'dih'+new_line
    defined_atoms.append(dih_start[3])
    #str(get_dihedral(rotamer.atom_list[dih_start[3]].get_xyz(), 
    #                 rotamer.atom_list[dih_start[2]].get_xyz(), 
    #                 rotamer.atom_list[dih_start[1]].get_xyz(), 
    #                 rotamer.atom_list[dih_start[0]].get_xyz()))
    dih_def = ''
    #if relax:
    dih_def='D 1 2 3 4 S '+str(set_steps)+' '+str(set_change)
    #else:
    #    dih_def='dih '+str(set_dih)+' '+str(set_steps)+' '+str(set_change)
    
    #now do all atoms attached to dih[0&3], and propigate out, then go back and do dih[1&2]attachments
    def add_zline(curr_num, trace):
        bond_num=trace[curr_num]
        ang_num=trace[bond_num]
        dih_num=trace[ang_num]
        add_string = z_name[curr_num]+' '+z_name[bond_num]+' '+str(get_dist(rotamer.atom_list[curr_num].get_xyz(), 
                                                              rotamer.atom_list[bond_num].get_xyz()))+' '
        add_string += z_name[ang_num]+' '+str(get_angle(rotamer.atom_list[curr_num].get_xyz(), 
                                                        rotamer.atom_list[bond_num].get_xyz(), 
                                                        rotamer.atom_list[ang_num].get_xyz()))+' '
        add_string += z_name[dih_num]+' '+str(get_dihedral(rotamer.atom_list[curr_num].get_xyz(), 
                                                        rotamer.atom_list[bond_num].get_xyz(), 
                                                        rotamer.atom_list[ang_num].get_xyz(), 
                                                        rotamer.atom_list[dih_num].get_xyz()))+new_line
        defined_atoms.append(curr_num)
        return add_string
    
    def trace_linker(to_add, trace, start_num):
        add_string = ''
        for num in to_add:
            trace[num] = start_num
            add_string+=add_zline(num, trace)
            pass_add = [x for x in rotamer.bond_dict[num] if x not in defined_atoms]
            if len(pass_add) >0:
                add_string+=trace_linker(pass_add, trace, num)
        return add_string
    
    #dih[0]
    traceback_dict = {}
    traceback_dict[dih_start[0]] = dih_start[1]
    traceback_dict[dih_start[1]] = dih_start[2]
    add_z = [x for x in rotamer.bond_dict[dih_start[0]] if x not in defined_atoms]
    out_string+=trace_linker(add_z, traceback_dict, dih_start[0])
    
    #dih[3]
    traceback_dict = {}
    traceback_dict[dih_start[3]] = dih_start[2]
    traceback_dict[dih_start[2]] = dih_start[1]
    add_z = [x for x in rotamer.bond_dict[dih_start[3]] if x not in defined_atoms]
    out_string+=trace_linker(add_z, traceback_dict, dih_start[3])
    
    #don't want to define atoms in the context of the bond that will be rotated if possible
    def check_for_dihedral(part_trace, fallback_trace, orig_atom):
        pos_start = [x for x in rotamer.bond_dict[part_trace[orig_atom]] if x != orig_atom]
        if len(pos_start) > 0:
            part_trace[part_trace[orig_atom]] = pos_start[0]
            return part_trace
        return fallback_trace
    
    #dih[1]
    fall_dict = {}
    #1->2->3
    fall_dict[dih_start[2]] = dih_start[3]
    fall_dict[dih_start[1]] = dih_start[2]
    part_dict={}
    part_dict[dih_start[1]] = dih_start[0]
    traceback_dict = check_for_dihedral(part_dict, fall_dict, dih_start[1])
    
    add_z = [x for x in rotamer.bond_dict[dih_start[1]] if x not in defined_atoms]
    out_string+=trace_linker(add_z, traceback_dict, dih_start[1])
    
    #dih[2]
    fall_dict = {}
    #2->1->0
    fall_dict[dih_start[2]] = dih_start[1]
    fall_dict[dih_start[1]] = dih_start[0]
    part_dict={}
    part_dict[dih_start[2]] = dih_start[3]
    traceback_dict = check_for_dihedral(part_dict, fall_dict, dih_start[2])
    
    add_z = [x for x in rotamer.bond_dict[dih_start[2]] if x not in defined_atoms]
    out_string+=trace_linker(add_z, traceback_dict, dih_start[2])
    
    
    out_string+='\n'+dih_def
    #print(out_string)
    #print(out_string.count('\n'))
    heading = ''
    #if relax:
    route = '%chk='+path+'.chk'+new_line+'%mem=8000MB'+new_line+'%nprocshared='+nproc+new_line
    route += '#'+theory+' scrf=smd nosymm opt=modredundant'
    title = 'Title Card Required'
    chrgspin = '0 1'
    heading = route + new_line + new_line + title + new_line + new_line + chrgspin + new_line
    #else:
    #    route = '%chk='+path+'.chk'+new_line+'%mem=8000MB'+new_line+'%nprocshared='+nproc+new_line
    #    route += theory+' scrf=smd nosymm scan'
    #    title = 'Title Card Required'
    #    chrgspin = '0 1'
    #    heading = route + new_line + new_line + title + new_line + new_line + chrgspin + new_line
    #print(heading+out_string)
    with open(path+'.gjf', 'w') as out_file:
        out_file.write(heading+out_string)

# Run main() from command line
if __name__ == '__main__':
    main ()