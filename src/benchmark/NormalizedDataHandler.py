import os

import pandas as pd

#class handling the normalized data

#dictionary from file name to content

#function to build classifiers for the datasets in dict


class NormalizedDataHandler:
    
    def __init__(self, file_list):
        data_dict = dict()
        for file in file_list:
            data_name = os.path.basename(file)
            data_content = pd.read_csv(file, sep = '\t')

        self.data = data_dict 
