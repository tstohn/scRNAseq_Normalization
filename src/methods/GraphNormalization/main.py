from NormalizationGraph import NormalizationGraph
import pandas as pd
import sys
import os.path

def read_data(data_file):
    data = pd.read_csv(data_file, sep='\t')
    return(data)

def normalize_data(dataset):
    G = NormalizationGraph(dataset)
    normalized_score_table = G.get_normalized_score()

    return(normalized_score_table)

def main():
    if(len(sys.argv) != 2):
        print("ERROR: use script with <python3 graphNormalization.py [directory of datasets]>\n")
        exit(1)
    data = read_data(sys.argv[1])
    result = normalize_data(data)

    result.to_csv("/Users/t.stohn/Desktop/Normalization/PIPELINE/scRNAseq_Normalization/bin/NORMALIZED_DATASETS/" + os.path.basename(os.path.dirname(sys.argv[1])) + "/GraphNormalized_x.tsv", sep='\t')

if __name__ == '__main__':
    main()
