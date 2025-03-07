from argparse import ArgumentParser

parser = ArgumentParser(description="command line arguments")
parser.add_argument('-i', '--input', required=True, help="Input file")
argument=parser.parse_args()

#filename = 'DHA_3D_Rotamer_top_scan_0.gjf'
filename = argument.input

out_str = ''
for i in range(1, 36):
    out_str = ''
    float_val = str(float(i)*10)
    int_val = str(i)
    with open(filename, 'r') as in_file:
        for line in in_file:
            line_split = line.split()
            if (len(line_split) > 0 and line_split[0] == 'C4') or line.strip() == 'D 5 1 2 3 0.0 F':
                temp = line.replace('0.0', float_val)
                out_str += temp
            elif '.chk' in line:
                temp = line.replace('0', int_val)
                out_str += temp
            else:
                out_str += line
    with open(filename.replace('0', int_val), 'w') as out_file:
        out_file.write(out_str)