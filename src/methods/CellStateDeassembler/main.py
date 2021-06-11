from NormalizationGraph import NormalizationGraph
import pandas as pd
import sys
import os.path
import re
from os import listdir, makedirs
from os.path import isfile, join

def read_data(data_file):
    data = pd.read_csv(data_file, sep='\t')
    return(data)

def ensure_dir(file_path):
    if not os.path.exists(file_path):
        makedirs(file_path)

def run_clusterAnalysis_on_dataset(data_dir, correlation, output_dir = ""):    
    data = read_data(data_dir)
    result =None

    if(output_dir):
        ensure_dir(output_dir)
        result.to_csv(output_dir + "/GraphNormalized.tsv", sep='\t')
    else:
        result.to_csv(os.path.dirname(data_dir) + "/GraphNormalized.tsv", sep='\t')

def main():
    if(len(sys.argv) != 3):
        print("ERROR: use script with <python3 CellStateDeassembler.py [directory of datasets] [output directory]>\n")
        exit(1)
    data_dir = sys.argv[1]
    output_dir = sys.argv[2]

    dataset_dir = "./bin/FILTERED_DATASETS"
    if(data_dir == "ALL"):
        #store all not ini files in list (store full path)
        dataset_list = list()
        for f in listdir(dataset_dir):
            f_path = join(dataset_dir, f)
            if isfile(f_path) and (f_path.endswith(".tsv")):
                dataset_list.append(f_path)

        #for every dataset in this list
        for dataset in dataset_list:
            output_dir = "bin/NORMALIZED_DATASETS/" + os.path.basename(os.path.splitext(dataset)[0])
            print("Running Graph normalization on: " + dataset)
            run_clusterAnalysis_on_dataset(dataset, output_dir)
    else:
        run_clusterAnalysis_on_dataset(data_dir)

if __name__ == '__main__':
    main()
