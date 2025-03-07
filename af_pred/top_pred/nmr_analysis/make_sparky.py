import argparse

aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
     'DAL': 'C', 'DBB':'C', 'ABA': 'C'} #conjugated resi more like cys

def parse_file(input_file, conn_file = ''):
    nmr_dict = {}
    loop = 0
    count = 0
    cys_list = []
    if conn_file != '':
        with open(conn_file, 'r') as in_file:
            for line in in_file:
                split = line.split(',')
                cys_list.append(split[0])
                cys_list.append(split[3])
    with open(input_file, 'r') as in_file:
        lines = in_file.readlines()
        length = len(lines)
        for line in lines:
            if 'loop_' in line:
                loop = count #sets number, not shared memory
                while '_Atom_chem_shift' in lines[loop +1]:
                    loop += 1
            count += 1 #will have index of current line
        #print(lines[loop], lines[loop+1])
        split_len = len(lines[loop+2].split()) #26
        i = loop+2
        #print(lines[i], lines[i+1])
        while len(lines[i].split()) == split_len:
            #print(lines[i])
            adj = 0
            if split_len == 24:
              adj = -1
            split = lines[i].split()
            index = split[6+adj]
            resi = split[7+adj]
            atom = split[8+adj]
            shift = split[11+adj]
            #print(index, resi, atom, shift)

            if index not in nmr_dict.keys():
                nmr_dict[index] = {}

                if index in cys_list:
                    nmr_dict[index]['AA'] = 'C'
                elif resi in aa.keys():
                    nmr_dict[index]['AA'] = aa[resi]
                else:
                    nmr_dict[index]['AA'] = 'A'

                nmr_dict[index]['HA'] = '0.00'
                nmr_dict[index]['H'] = '0.00'
                nmr_dict[index]['CA'] = '0.00'
                nmr_dict[index]['C'] = '0.00' #verify this if possible
                nmr_dict[index]['CB'] = '0.00'
                nmr_dict[index]['N'] = '0.00'

            #if index in cys_list:
            #    nmr_dict[index][atom] = shift
            if resi in aa.keys():
                nmr_dict[index][atom] = shift

            i += 1
    #print nmr_dict
    return nmr_dict

def sparky(filename, spark_dict):
    out_name = filename.replace('.txt', '_spark.txt')
    if filename != out_name:
        with open(out_name, 'w') as out_file:
            out_file.write('#NUM AA HA HN N15 CA CB CO\n')
            keys = spark_dict.keys()
            keys = [int(x) for x in keys]
            keys.sort()
            keys = [str(x) for x in keys]
            for key in keys:
                out_str = key + ' '
                #print(key, spark_dict[key])
                for item in ['AA', 'HA', 'H', 'N', 'CA', 'CB', 'C']:
                    out_str += spark_dict[key][item] + ' '
                out_file.write(out_str[:-1]+'\n')
                

def args():
    parser = argparse.ArgumentParser(description="Process NMR data")
    parser.add_argument("-i", "--input", required=True, help="Input file containing NMR data")
    parser.add_argument("-c", "--connect", default='', help="Optional file specifying lanthipeptide connectivity")
    
    args = parser.parse_args()
    return args.input, args.connect
    
    #parse_file(input_file)

if __name__ == "__main__":
    in_file, conn_file = args()
    data_dict = parse_file(in_file, conn_file)
    sparky(in_file, data_dict)
 
