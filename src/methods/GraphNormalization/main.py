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

def parse_ini_for_batch_effect_removal(data_dir):
    #check path to ini files for datasets from settings.ini
    settings = open("settings.ini", "r")
    path = ""
    while settings:
        line = settings.readline()
        print(line)
        data_path = re.match(("datasets=(.*)"), line)
        if(data_path):
            path = data_path[1].rstrip("\n")
            break
    #parse ini file for dataset
    data_ini = open(path + "/" + os.path.basename(os.path.splitext(data_dir)[0]) + ".ini", "r")
    while data_ini:
        line = data_ini.readline()
        batch_eff = re.match(("removeBatchEffect=(\\d)"), line)
        if(batch_eff):
            return(int(batch_eff[1].rstrip("\n")))

def normalize_data(dataset, correlation, remove_batch_effect=1):
    corr_method = "spearman"
    corr_threshold=correlation
    G = NormalizationGraph(dataset, corr_method, corr_threshold)
    p_val=0.05
    cohend_val=0.5
    take_log=False
    if(remove_batch_effect):
        G.remove_batch_effect()
    G.normalize_by_library_size()
    G.get_normalized_score(p_val, cohend_val, take_log)
    return(G.data)

def ensure_dir(file_path):
    if not os.path.exists(file_path):
        makedirs(file_path)

def run_normalization_on_dataset(data_dir, correlation, output_dir = ""):    
    data = read_data(data_dir)
    remove_batch_effect = parse_ini_for_batch_effect_removal(data_dir)
    try:
        result = normalize_data(data, correlation, remove_batch_effect)
    except:
        print("Could not normalize dataset " + data_dir)
        return
    if(output_dir):
        ensure_dir(output_dir)
        result.to_csv(output_dir + "/GraphNormalized.tsv", sep='\t')
    else:
        result.to_csv(os.path.dirname(data_dir) + "/GraphNormalized.tsv", sep='\t')

def main():
    if(len(sys.argv) != 3 and len(sys.argv) != 2):
        print("ERROR: use script with <python3 graphNormalization.py [directory of datasets]>\n")
        exit(1)
    data_dir = sys.argv[1]
    if(len(sys.argv) == 3):
        correlation = float(sys.argv[2])
    else:
        correlation = 0.7
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
            run_normalization_on_dataset(dataset, correlation, output_dir)
    else:
        run_normalization_on_dataset(data_dir, correlation)

if __name__ == '__main__':
    main()
