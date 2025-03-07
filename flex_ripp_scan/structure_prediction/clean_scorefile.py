import argparse

def clean_scorefile(input_file, output_file):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Ensure there are at least two lines for SEQUENCE and total_score
    if len(lines) < 2:
        print("Error: Scorefile does not contain enough lines.")
        return

    # Keep the first line (SEQUENCE) and second line (total_score)
    sequence_line = lines[0]
    total_score_line = lines[1]
    num_elements = len(total_score_line.split())

    # Write to the output file
    with open(output_file, 'w') as outfile:
        outfile.write(sequence_line)
        outfile.write(total_score_line)

        # Keep any other lines that have the same number of elements as the total_score line
        for line in lines[2:]:
            if len(line.split()) == num_elements:
                outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description='Clean a Rosetta scorefile.')
    parser.add_argument('-s', '--scorefile', required=True, help='Path to the scorefile to be cleaned.')

    args = parser.parse_args()

    # Construct the output filename by adding '_clean.sc' to the base name
    base_name = args.scorefile.split('.')[0]
    output_file = f"{base_name}_clean.sc"

    # Clean the scorefile
    clean_scorefile(args.scorefile, output_file)
    print(f"Cleaned scorefile saved as: {output_file}")

if __name__ == '__main__':
    main()
