import argparse
import pandas as pd

def process_scorefile(scorefile):
    # Load the scorefile into a pandas DataFrame, assuming space-separated values
    df = pd.read_csv(scorefile, delim_whitespace=True, skiprows=1)

    # Check required columns
    required_columns = ['total_score', 'rmsd_bb', 'description', 'nmr_cst']
    if not all(col in df.columns for col in required_columns):
        raise ValueError(f"Scorefile must contain the columns: {', '.join(required_columns)}")

    # Extract the base descriptions (everything before "_xxxx")
    df['base'] = df['description'].str.extract(r'(.*)_\d{4}')

    # Group by base and find the row with the lowest total_score for each base
    result = df.loc[df.groupby('base')['nmr_cst'].idxmin(), ['rmsd_bb', 'total_score', 'nmr_cst', 'description']]

    return result

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(description="Process a Rosetta scorefile to identify structures that coorespond with low nmr_cst violations, report good scoring ones.")
    parser.add_argument('-i', '--in_file', help='Input scorefile', required=True)

    # Parse arguments
    args = parser.parse_args()

    base = args.in_file.split('/')[-1].split('_')[0]

    # Set pandas display options
    pd.set_option('display.max_columns', None)  # Show all columns
    pd.set_option('display.max_rows', None)     # Show all rows
    pd.set_option('display.max_colwidth', None) # Show full width of each column
    pd.set_option('display.width', None)        # Auto-detect the display width

    try:
        # Process the scorefile
        result = process_scorefile(args.in_file)

        #for i in range(0, 11):
        #    quantile = float(i)/10
        #    nmr_cst_thresh = round(result['nmr_cst'].quantile(quantile), 2)
        #    temp_df = result[result['nmr_cst'] <= nmr_cst_thresh ]
        #    print(quantile, nmr_cst_thresh, round(temp_df['total_score'].quantile(0.0), 2))

        df = result.nsmallest(n=200, columns=['nmr_cst'])
        df = df.nsmallest(n=20, columns=['total_score'])
        #print(df)
        with open(f'{base}_nmr_opt.txt', 'w') as out_file:
            for des in df['description']:
                out_file.write(f'{des}.pdb.gz\n')

        # Save the result to the output file
        #out_name = args.in_file.replace('.sc', '') + '.dat'
        #result.to_csv(out_name, sep=' ', index=False, header=True)
        #print(f"Results saved to {out_name}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
