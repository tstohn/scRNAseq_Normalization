
# Import libraries
from scvi.dataset import CsvDataset
import os
import numpy as np
import pandas as pd
from scvi.models import *
from scvi.inference import UnsupervisedTrainer
from scvi.dataset.dataset import GeneExpressionDataset
import torch
import sys
import os.path
from os import listdir, makedirs
from os.path import isfile, join

def add_scVInormalized_data(data, normData):
    melted_normData = pd.melt(normData, id_vars=['sample_id'], var_name='ab_id', value_name='ab_count_normalized')
    merged_data = data.merge(melted_normData, on=['sample_id', 'ab_id'], how='left')
    return(merged_data)

def main():
    #params
    if(len(sys.argv) != 4):
        print("ERROR: use script with <python3 SanityNormalization/main.py [directory of dataset] [output directory]>\n")
        exit(1)
    data_dir = sys.argv[1]
    output_dir = sys.argv[2]

    #read in data
    dataFrame = pd.read_table(data_dir)
    selected_columns = dataFrame[['sample_id', 'ab_id', 'ab_count']]
    pivoted_df = selected_columns.pivot(index='sample_id', columns='ab_id', values='ab_count')
    pivoted_df.reset_index(inplace=True)
    indices = pivoted_df["sample_id"]
    abs = pivoted_df.columns
    pivoted_df.drop(['sample_id'], axis=1, inplace=True)

    #perform normalization

    #WRITE RESULTS
    ensure_dir(output_dir)
    result.to_csv(output_dir + "/scVI.tsv", sep='\t')

if __name__ == '__main__':
    main()
