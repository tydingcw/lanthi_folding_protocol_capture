import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description = 'take in a score file and output the lowest scoring ddg structures')
parser.add_argument('-i', '--in_file', help = 'input scorefile', required = True)
#parser.add_argument('-n', '--number', help = 'number of struct to take', default=10)
#parser.add_argument('-c', '--column', help = 'column of interest', default='ddg_des')
#parser.add_argument('-o', '--out', help = 'output type: print df or desc default is df', default='df')
#parser.add_argument('-s', '--scorefile', help = 'write out scorefile info for clustering', action='store_true', default=False)
parser.add_argument('-c', '--cutoff', help = 'score cutoff', type=float, default=50.0)
argument=parser.parse_args()

#filename = sys.argv[1]
#argument.number = int(argument.number)

cols = None
with open(argument.in_file, 'r') as in_file:
    first_line = in_file.readline()
    first_line = in_file.readline()
    cols=first_line.split()

basename = argument.in_file.split('.')[0]

#print(cols)

#if argument.scorefile:
#    argument.out = 'df'


df = pd.read_csv(argument.in_file, delim_whitespace=True, skiprows=1)[cols]
df = df[df['total_score'] < argument.cutoff]
#print(df.columns)
data_cols = [x for x in cols if 'ring' in x]
data_cols.append('rmsd_bb')
data_cols.append('rmsd_full')
data_cols.append('nmr_cst')
data_cols.append('SASA')
eng_cols = [x for x in data_cols if 'eng' in x]
for col in eng_cols:
    data_cols.remove(col)
    
#print(data_cols)
#if argument.out == 'df':
for col in data_cols:
    eng_col = 'total_score'
    if 'ring' in col:
        eng_col =  col.split('_')[0] + '_eng'
    out_df = df[[col, eng_col]]
    out_path = basename + '_' + col + '.dat'
    out_df.to_csv(out_path, sep=' ', index=False)
    col_upper = col.upper()
    #out_df = out_df[df[argument.score] < argument.cutoff]
    #print(f'bash pnear.sh {out_path} {eng_col} {col}')
    os.system(f'bash pnear.sh {out_path} {eng_col} {col} {col_upper}')
    
    #df = df.nsmallest(n=argument.number, columns=[argument.column])
    #if argument.scorefile:
        #out_str = ''
        #file_list = ''
        #for index, row in df.iterrows():
        #    #print(row['description'].split('/')[-1])
        #    out_str += row['description'].split('/')[-1]+'.pdb.gz ' + str(row[argument.column]) + '\n'
        #    file_list += row['description']+'.pdb.gz '
        #with open(argument.in_file.split('.')[0] + '_extract.txt', 'w') as out_file:
        #    out_file.write(out_str)
        #with open(argument.in_file.split('.')[0] + '_extract.options', 'w') as out_file:
        #    out_file.write('-in:file:s ' + file_list + '\n')
        #    out_file.write(options)
        #    out_file.write(argument.in_file.split('.')[0] + '_extract.txt')
        #print(out_str)
    #else:
    #    print(df[['description', argument.column]])
#elif argument.out == 'desc':
#    for des in df.nsmallest(n=argument.number, columns=[argument.column])['description']:
#        print(des)
