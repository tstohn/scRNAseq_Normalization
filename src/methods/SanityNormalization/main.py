
# Import libraries
import os
import numpy as np
import pandas as pd
import torch
import sys
import os.path
from os import listdir, makedirs
from os.path import isfile, join
import subprocess
import shutil
import re

def add_Sanitynormalized_data(data, normData):
    #rename ab column (renamed to run Sanity)
    normData.rename(columns={"GeneID": "ab_id"}, inplace = True)
    #melt into tibble format
    melted_normData = pd.melt(normData, id_vars=['ab_id'], var_name='sample_id', value_name='ab_count_normalized')
    #merge with origional data
    merged_data = data.merge(melted_normData, on=['sample_id', 'ab_id'], how='left')
    return(merged_data)

def main():
    #params // ensure directories exist
    if(len(sys.argv) != 3):
        print("ERROR: use script with <python3 SanityNormalization/main.py [input file] [output directory]>\n")
        exit(1)
    data_file = sys.argv[1] #file which should be normalized (tibble format)
    output_dir = sys.argv[2] #directory where all final normalized tables should be stored
    if not os.path.exists(output_dir):
        makedirs(output_dir)  
    #folder for Sanity tmp results
    sanityTmp = "./src/methods/SanityNormalization/Sanity/Tmp"
    if not os.path.exists(sanityTmp):
        os.makedirs(sanityTmp) 

    #read in data
    dataFrame = pd.read_table(data_file)
    selected_columns = dataFrame[['sample_id', 'ab_id', 'ab_count']]
    pivoted_df = selected_columns.pivot(index='ab_id', columns='sample_id', values='ab_count')
    pivoted_df.reset_index(inplace=True)
    ab_ids = pivoted_df["ab_id"]
    sample_ids = pivoted_df.columns
    pivoted_df.rename(columns={"ab_id": "GeneID"}, inplace = True)

    #perform normalization
    #write data to tmp file, so that Sanity can handle it
    pivoted_df.to_csv(sanityTmp + "/SanityTmpMatrix.tsv", sep='\t', index=False)
    commandString = "./src/methods/SanityNormalization/Sanity/bin/Sanity -f " + sanityTmp + "/SanityTmpMatrix.tsv -d " + sanityTmp + "" #additional parameters
    try:
        subprocess.run([commandString], shell = True, check = True)
    except Exception as e: 
        print(e)
                  
    #read in tmp solution of sanity // add ab_count_normalized to data frame
    tmpSanitySolution = pd.read_table(sanityTmp + "/log_transcription_quotients.txt")
    result = add_Sanitynormalized_data(dataFrame, tmpSanitySolution)

    #WRITE RESULTS 
    result.to_csv(output_dir + "/Sanity.tsv", sep='\t')

    #clean up (tmp files with Sanity results, pivoted table)
    shutil.rmtree(sanityTmp)

if __name__ == '__main__':
    main()
