#Clay Tydings
import argparse
import pandas as pd

def get_ensemble_list(filename):
    ens_structs = []
    with open(filename, 'r') as in_file:
        for raw_line in in_file:
            line = raw_line.strip()
            if 'score' not in line:
                ens_structs.append(line.split('/')[-1].split('.')[0])
    return ens_structs

def get_ens_min_scoring(scorefile, ensemble):
    cols = None
    with open(scorefile, 'r') as in_file:
        first_line = in_file.readline()
        first_line = in_file.readline()
        cols=first_line.split()
    
    full_df = pd.read_csv(scorefile, delim_whitespace=True, skiprows=1)[cols]
    
    scores = {}
    
    #print(ensemble)
    #for name in ensemble:
    #    df = full_df.loc[full_df['description'] == name ]
    for index, row in full_df.iterrows():
        #print()
        if row['description'].split('/')[-1].split('.')[0] in ensemble:
            scores[float(row['total_score'])] = row['description']
            
    #return df.nsmallest(n=1, columns=['total_score'])['description'].iloc[0]
    low_score = min (scores.keys())
    return scores[low_score]

def main():
    parser = argparse.ArgumentParser(description='get the lowest scoring member of the ensemble')
    parser.add_argument('-s', '--scorefile', help='scorefile')
    parser.add_argument('-e', '--ensemble', help='file with ensemble members')
    
    args = parser.parse_args()
    
    print(get_ens_min_scoring(args.scorefile, get_ensemble_list(args.ensemble)))


if __name__ == '__main__':
    main()
