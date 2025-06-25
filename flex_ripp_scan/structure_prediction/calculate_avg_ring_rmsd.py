#Clay Tydings
import argparse
import pandas as pd

#def get_ensemble_list(filename):
#    ens_structs = []
#    with open(filename, 'r') as in_file:
#        for raw_line in in_file:
#            line = raw_line.strip()
#            if 'score' not in line:
#                ens_structs.append(line.split('/')[-1].split('.')[0])
#    return ens_structs

def get_ring_avg(scorefile):
    cols = None
    with open(scorefile, 'r') as in_file:
        first_line = in_file.readline()
        first_line = in_file.readline()
        cols=first_line.split()
    
    full_df = pd.read_csv(scorefile, delim_whitespace=True, skiprows=1)[cols]
    
    scores = {}
    
    eng_cols = [x for x in cols if 'eng' in x and 'ring' in x]
    rmsd_cols = [x for x in cols if 'rmsd_bb' in x and 'ring' in x]

    #to record the eng and rmsd data
    data = []
    data.append(['avg_rmsd', 'avg_energy'])

    for index, row in full_df.iterrows():
        eng_val = 0
        rmsd_val = 0
        for col in eng_cols:
            eng_val += float(row[col])
        for col in rmsd_cols:
            rmsd_val += float(row[col])
        eng_val /= len(eng_cols)
        rmsd_val /= len(rmsd_cols)
        if pd.notna(eng_val) and pd.notna(rmsd_val):
            data.append([str(rmsd_val), str(eng_val)])
 
    return data

def main():
    parser = argparse.ArgumentParser(description='generate dat file for ring averaged rmsd and energy')
    parser.add_argument('-s', '--scorefile', help='scorefile')
    
    args = parser.parse_args()
    
    rmsd_data = get_ring_avg(args.scorefile)
    outname = args.scorefile.replace('.sc', '') + '_ring_avg_rmsd.dat'
    with open(outname, 'w') as csvfile:
        for row in rmsd_data:
            csvfile.write(' '.join(row) + '\n')


if __name__ == '__main__':
    main()
