from NormalizationGraph import NormalizationGraph
import pandas as pd
import sys
import os.path
from os import listdir, makedirs
from os.path import isfile, join

def read_data(data_file):
    data = pd.read_csv(data_file, sep='\t')
    filename = os.path.splitext(data_file)[0] + ".ini"
    col_name_dict = {}
    with open(filename, 'r') as f:
        line = f.readline()
        while (line):
            if (line.startswith("sample_id")):
                split=line.split("=")
                from_val = split[1].rstrip()
                to_val = split[0]
                col_name_dict[from_val]=to_val
            elif(line.startswith("ab_id")):
                split=line.split("=")
                from_val = split[1].rstrip()
                to_val = split[0]
                col_name_dict[from_val]=to_val
            elif(line.startswith("ab_count")):
                split=line.split("=")
                from_val = split[1].rstrip()
                to_val = split[0]
                col_name_dict[from_val]=to_val
            elif(line.startswith("batch_id")):
                split=line.split("=")
                from_val = split[1].rstrip()
                to_val = split[0]
                col_name_dict[from_val]=to_val                
            elif(line.startswith("cluster_id")):
                split=line.split("=")
                from_val = split[1].rstrip()
                to_val = split[0]
                col_name_dict[from_val]=to_val
            elif(line.startswith("ab_type")):
                split=line.split("=")
                from_val = split[1].rstrip()
                to_val = split[0]
                col_name_dict[from_val]=to_val
            line = f.readline()
    for key in col_name_dict:
        print(key)
        data.rename(columns = {key: col_name_dict[key]}, inplace = True)

    print(data)
    return(data)

def normalize_data(dataset):
    G = NormalizationGraph(dataset)
    normalized_score_table = G.get_normalized_score()
    return(normalized_score_table)

def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

def run_normalization_on_dataset(data_dir, output_dir = ""):    
    data = read_data(data_dir)
    result = normalize_data(data)
    if(output_dir):
        ensure_dir(output_dir)
        result.to_csv(output_dir + "/GraphNormalized.tsv", sep='\t')
    else:
        result.to_csv(os.path.dirname(data_dir) + "/GraphNormalized.tsv", sep='\t')

def main():
    if(len(sys.argv) != 2):
        print("ERROR: use script with <python3 graphNormalization.py [directory of datasets]>\n")
        exit(1)
    data_dir = sys.argv[1]
    dataset_dir = ""
    if(data_dir == "ALL"):
        #parse directory name
        filename = "settings.ini"
        with open(filename, 'r') as f:
            while (1):
                line = f.readline()
                if (line.startswith("datasets")):
                    dataset_dir = line.split("=")[1]
                    dataset_dir=dataset_dir.rstrip()
                    break
        assert(dataset_dir != "")

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
            run_normalization_on_dataset(dataset, output_dir)
    else:
        run_normalization_on_dataset(data_dir)

if __name__ == '__main__':
    main()
