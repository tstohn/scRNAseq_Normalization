import sys
import os
from os import listdir
from os.path import isfile, join, isdir

def parse_error(dir):
    print("ERROR in " + dir + ": directory must contain only noralized files or directories of normalized files")
    exit(1)

def listdir_nohidden(dir):
    for file in os.listdir(dir):
        if not (file.startswith('.') or file.startswith('_')):
            yield file

def load_datasets_for_benchmark(dir):
    datasets = list()

    #get all elements in directory
    dir_content = listdir_nohidden(dir)
    one_dataset_only = True
    file_list = list()
    for sub_dir in dir_content:
        sub_dir_path = join(dir,sub_dir)
        if(os.path.isdir(sub_dir_path) and not file_list):
            one_dataset_only = False
            for file in listdir_nohidden(sub_dir_path):
                if(isfile(join(sub_dir_path, file)) and not os.path.basename(file).endswith("metaData.tsv")):
                    file_list.append(join(sub_dir_path, file))
            datasets.append(file_list.copy())
            file_list.clear()
        elif(one_dataset_only):
            file_path = join(dir, sub_dir)
            if(isfile(file_path) and not os.path.basename(file_path).endswith("metaData.tsv")):
                file_list.append(file_path)
        else:
            parse_error(dir)
    if(file_list):
        datasets.append(file_list.copy())


    print(datasets)
    return datasets