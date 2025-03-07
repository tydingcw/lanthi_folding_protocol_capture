from argparse import ArgumentParser
import numpy as np

class Rotamer:
    def __init__(self):
        self.atom_list = []
        self.bond_dict = {}
        self.sulfur_pos = -1
        self.conf_score = None
        self.clash_score = None
    
    def add_atom(self, atom):
        self.atom_list.append(atom)
        if atom.name == 'S':
            self.sulfur_pos = len(self.atom_list)-1
    
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

def get_dist(atom1, atom2):
    dist = ((atom1[0]-atom2[0])**2 + (atom1[1]-atom2[1])**2 + 
            (atom1[2]-atom2[2])**2)**(0.5)
    return dist

parser = ArgumentParser(description="command line arguments")
parser.add_argument('-i', '--input', required=True, help="Input file")
parser.add_argument('-m', '--min', default='min', help="Input file")
args=parser.parse_args()

def get_min(filename):
    #rotamer_list = []
    out_str = ''
    temp_str = ''
    #read_frame = False
    get_conf_score = False
    conf_score = 1000.0
    best_conf = False
    #get_clash_score = False
    with open(filename, 'r') as file:
        for raw_line in file:
            #line = raw_line.strip()
            #line_list = line.split()
            #if len(line_list)==11 and not read_frame:
            #    read_frame = True
            #    rotamer_list.append(Rotamer())
            #elif 'END' in raw_line:
            #    read_frame = False
            #elif read_frame and len(line_list) == 16:
            #    #have an atom
            #    new_atom = Atom(line_list[3], line_list[0], 
            #                   line_list[1], line_list[2])
            #    rotamer_list[len(rotamer_list)-1].add_atom(new_atom)
            #elif read_frame and len(line_list) == 7:
            #    #have a bond
            #    rotamer_list[len(rotamer_list)-1].add_bond(line_list[0], 
            #                                               line_list[1])
            temp_str += raw_line
            if get_conf_score:
                get_conf_score = False
                if float(raw_line) < conf_score:
                    conf_score = float(raw_line)
                    best_conf = True
            elif '<ConfScore>' in raw_line:
                get_conf_score = True
            elif '$$$$' in raw_line:
                if best_conf:
                    out_str = temp_str[:-1]
                    best_conf = False
                temp_str = ''
            #if get_clash_score:
            #    get_clash_score = False
            #    rotamer_list[len(rotamer_list)-1].clash_score = float(line)
            #elif '<ConfClashScore>' in line:
            #    get_clash_score = True
            
        #return rotamer_list
    file_split = filename.split('.')
    with open(file_split[0]+'_top.sdf', 'w') as out_file:
        out_file.write(out_str)

def get_bin(rotamer):
    #phi C-N-CA-C (18,14,6,13), psi N-CA-C-N (14,6,13,16)
    phi = get_dihedral(rotamer.atom_list[17].get_xyz(), rotamer.atom_list[13].get_xyz(), rotamer.atom_list[5].get_xyz(), rotamer.atom_list[12].get_xyz())
    psi = get_dihedral(rotamer.atom_list[13].get_xyz(), rotamer.atom_list[5].get_xyz(), rotamer.atom_list[12].get_xyz(), rotamer.atom_list[15].get_xyz())
    if phi > 0 and psi > 0:
        #print(phi, psi, 'X')
        return 'X' #mirror of A
    elif phi <= 0 and psi <= 0:
        #print(phi, psi, 'A')
        return 'A' #alpha helix
    elif phi <= 0 and psi > 0:
        #print(phi, psi, 'B')
        return 'B' #beta sheet
    elif phi > 0 and psi <= 0:
        #print(phi, psi, 'Y')
        return 'Y' #mirror of B

def get_dmin(filename):
    rotamer_list = [] 
    #list_alpha = [] 
    #list_beta = [] 
    out_str_alpha = '' #for d amino acid, phi is pos, psi is pos
    out_str_beta = ''#for d amino acid, phi is pos, psi is neg
    temp_str = ''
    read_frame = False
    get_conf_score = False
    conf_alpha_score = 1000.0
    conf_beta_score = 1000.0
    #best_conf = False
    best_alpha_conf=False
    best_beta_conf=False
    #get_clash_score = False
    with open(filename, 'r') as file:
        for raw_line in file:
            temp_str += raw_line
            line = raw_line.strip()
            line_list = line.split()
            if len(line_list)==11 and not read_frame:
                read_frame = True
                rotamer_list.append(Rotamer())
            elif 'END' in raw_line:
                read_frame = False
            elif read_frame and len(line_list) == 16:
                #have an atom
                new_atom = Atom(line_list[3], line_list[0], 
                               line_list[1], line_list[2])
                rotamer_list[-1].add_atom(new_atom)
            elif read_frame and len(line_list) == 7:
                #have a bond
                rotamer_list[-1].add_bond(line_list[0], line_list[1])
                
            #rotamer will be complete once hit <ConfScore>, $$$$ is the end of the conf            
            if get_conf_score:
                get_conf_score = False
                if get_bin(rotamer_list[-1]) == 'X':
                    if float(raw_line) < conf_alpha_score:
                        conf_alpha_score = float(raw_line)
                        best_alpha_conf = True
                elif get_bin(rotamer_list[-1]) == 'Y':
                    if float(raw_line) < conf_beta_score:
                        conf_beta_score = float(raw_line)
                        best_beta_conf = True
            elif '<ConfScore>' in raw_line:
                get_conf_score = True
            elif '$$$$' in raw_line:
                if best_alpha_conf:
                    out_str_alpha = temp_str[:-1]
                    best_alpha_conf = False
                elif best_beta_conf:
                    out_str_beta = temp_str[:-1]
                    best_beta_conf = False
                temp_str = ''
            #if get_clash_score:
            #    get_clash_score = False
            #    rotamer_list[len(rotamer_list)-1].clash_score = float(line)
            #elif '<ConfClashScore>' in line:
            #    get_clash_score = True
            
        #return rotamer_list
    file_split = filename.split('.')
    with open(file_split[0]+'_top_dalpha.sdf', 'w') as out_file:
        out_file.write(out_str_alpha)
    with open(file_split[0]+'_top_dbeta.sdf', 'w') as out_file:
        out_file.write(out_str_beta)

if args.min == 'min': 
    get_min(args.input)
elif args.min == 'dmin': 
    print('not debugged yet')
    get_dmin(args.input)