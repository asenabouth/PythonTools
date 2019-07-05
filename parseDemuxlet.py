# Import libaries
import argparse
import pandas as pd
from glob import glob

def parse_individuals(x, split_indices = tuple()):
    # Store individuals
    individuals = []
    individual_boundaries = []
    n_individuals = len(split_indices)

    # Set boundaries
    start = 0
    end = len(x)

    for index in split_indices:
        # Get position of the index
        index_pos = split_indices.index(index) + 1
        
        # If it is less than the end of the list, set it as the end
        if (index_pos < n_individuals):
            individual_set = (start, index + 1)
            individual_boundaries.append(individual_set)
            start = index + 1
        else:
            individual_set = (start, end)
            individual_boundaries.append(individual_set)

    for start, end in individual_boundaries:
        individual = x[start:end]
        individual = ("-").join(individual)
        individuals.append(individual)
    return(individuals)

def parse_string(x):
    # Extract droplet type and strip from name string
    droplet_type = x.split("-")[0]
    x = x.replace((droplet_type + "-"), "")
    
    if (droplet_type in ["DBL", "AMB"]):
        x = x.split("-")
        split_indices = [n for n in range(0, len(x)) if ".CEL" in x[n]]
        individuals = parse_individuals(x, split_indices = split_indices)
        individuals = (",").join(individuals)
    else:
        individuals = x
    return(droplet_type, individuals)

def parseArguments():
    # Parser for inputs
    parser = argparse.ArgumentParser(description = "Extract assignments from 'demuxlet' best output.")
    parser.add_argument('-i', '--input', type = str, help = ".best file output from demuxlet.", required = True)
    parser.add_argument('-o', '--output', type = str, help = "Filename of output tsv file.", required = True)
    args = parser.parse_args()
    return(args.input, args.output)

if __name__ == "__main__":
    # Read in demuxlet output
    input_path, output_path = parseArguments()
    result_df = pd.read_csv(input_path, sep = "\t")

    # Parse best
    best = result_df["BEST"]
    parsed_strings = [parse_string(string) for string in best]
    droplet_types = [string[0] for string in parsed_strings]
    predictions = [string[1] for string in parsed_strings]

    # Add to data frame
    result_df["DROPLET.TYPE"] = droplet_types
    result_df["INDIVIDUALS"] = predictions

    output_column = ['BARCODE', 'RD.TOTL', 'RD.PASS', 'RD.UNIQ', 'N.SNP', 'DROPLET.TYPE', 'INDIVIDUALS', 'BEST', 'SNG.1ST',
        'SNG.LLK1', 'SNG.2ND', 'SNG.LLK2', 'SNG.LLK0', 'DBL.1ST', 'DBL.2ND',
        'ALPHA', 'LLK12', 'LLK1', 'LLK2', 'LLK10', 'LLK20', 'LLK00', 'PRB.DBL',
        'PRB.SNG1']

    result_df = result_df[output_column]
    result_df.to_csv(output_path, index = False, sep = "\t")
    print("%s parsed." % input_path)
        

