import os
import argparse

parser = argparse.ArgumentParser(description = 'take in a pdb code and find prepare the options files for generating 200,000 structures')
parser.add_argument('-i', '--input', help = 'input pdb code', required = True)
argument=parser.parse_args()

nstruct = len(os.listdir(argument.input)) #6VLJ
print(f'{argument.input} has {nstruct} structs')

multi_factor = int((200000-nstruct)/nstruct) #+ 1
num_struct = multi_factor * 100
#print(int((200000-nstruct)/nstruct)) #+ 1
#print(multi_factor, num_struct)
#print(multi_factor * 100 * 18)

if multi_factor > 1:
    options = [x for x in os.listdir('options/') if argument.input in x]
    for opt in options:
        temp_str = ''
        with open(f'options/{opt}', 'r') as in_option:
            for raw_line in in_option:
                temp_str += raw_line
        temp_str = temp_str.replace('100', str(num_struct)).replace('init', 'run1')
        with open(f'options_refined/{opt}', 'w') as out_option:
            out_option.write(temp_str)

#successful_runs = os.listdir(argument.input) #6VLJ
#run_dict = {}
#for filename in successful_runs:
#    num = filename.split('_')[0]
#    if num not in run_dict.keys():
#        run_dict[num] = 1
#    else:
#        run_dict[num] += 1
#
#good_list = []
#for key in run_dict.keys():
#    if run_dict[key] >= 19:
#        good_list.append(key)
#num_struct = int(50000/len(good_list))
#
#for item in good_list:
#    temp_str = ''
#    with open(f'options/{argument.input}_{item}.options', 'r') as in_option:
#        for raw_line in in_option:
#            temp_str += raw_line
#    temp_str = temp_str.replace('20', str(num_struct)).replace('init', 'run1')
#    with open(f'options_refined/{argument.input}_{item}.options', 'w') as out_option:
#        out_option.write(temp_str)
