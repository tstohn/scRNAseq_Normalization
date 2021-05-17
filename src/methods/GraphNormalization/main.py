from NormalizationGraph import NormalizationGraph
import pandas as pd
import sys

def read_data(data_file):
    data = pd.read_csv(data_file, sep='\t')
    return(data)

def normalize_data(dataset):
    G = NormalizationGraph(dataset)
    normalized_scores = G.get_normalized_score()



def main():
    if(len(sys.argv) != 2):
        print("ERROR: use script with <python3 graphNormalization.py [directory of datasets]>\n")
        exit(1)
    data = read_data(sys.argv[1])
    normalize_data(data)

if __name__ == '__main__':
    main()
