
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
    #merge with origional data (make surecolumn are of same data type)
    data.sample_id = data.sample_id.astype(str)
    melted_normData.sample_id = melted_normData.sample_id.astype(str)
    data.ab_id = data.ab_id.astype(str)
    melted_normData.ab_id = melted_normData.ab_id.astype(str)
    merged_data = data.merge(melted_normData, on=['sample_id', 'ab_id'], how='left')
    return(merged_data)

def main():

    if(len(sys.argv) != 3):
        print("ERROR: use script with <python3 scVINormalization/main.py [directory of dataset] [output directory]>\n")
        exit(1)
    data_dir = sys.argv[1]
    output_dir = sys.argv[2]
    if not os.path.exists(output_dir):
        makedirs(output_dir) 

    #read in data
    dataFrame = pd.read_table(data_dir)
    selected_columns = dataFrame[['sample_id', 'ab_id', 'ab_count']]
    pivoted_df = selected_columns.pivot(index='sample_id', columns='ab_id', values='ab_count')
    pivoted_df.reset_index(inplace=True)

    indices = pivoted_df["sample_id"]
    pivoted_df.drop(['sample_id'], axis=1, inplace=True)
    abs = pivoted_df.columns

    #write into GeneExpressionDataset (cell*gene matrix)
    dataset = GeneExpressionDataset()
    dataset.populate_from_data(
        X=pivoted_df.values,  # Use the values of the DataFrame as the expression matrix
        batch_indices=None,  # If you have batch information, you can specify it here
        gene_names=pivoted_df.columns,  # Gene names from the DataFrame columns
        cell_types=None  # Cell type information (optional)
    )

    # Process data (from https://github.com/YosefLab/scVI/blob/master/tests/notebooks/basic_tutorial.ipynb)
    use_batches=False
    use_cuda=True
    n_epochs=400
    vae = VAE(dataset.nb_genes, n_batch=dataset.n_batches * use_batches)
    trainer = UnsupervisedTrainer(vae, dataset, train_size=0.75, use_cuda=use_cuda)
    trainer.train(n_epochs=n_epochs)
    full = trainer.create_posterior(trainer.model, dataset, indices=np.arange(len(dataset)))
    #imputed_values = full.sequential().imputation()
    normalized_values = full.sequential().get_sample_scale()

    normalizedFrame = pd.DataFrame(normalized_values, columns = abs)
    normalizedFrame["sample_id"] = indices

    # Write output matrix
    result = add_scVInormalized_data(dataFrame, normalizedFrame)
    result.to_csv(output_dir + "/scVI.tsv", sep='\t')

if __name__ == '__main__':
    main()
