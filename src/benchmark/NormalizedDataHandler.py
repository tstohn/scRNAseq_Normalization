import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from random import randint

from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score, KFold, GridSearchCV, RandomizedSearchCV
from sklearn.tree import DecisionTreeClassifier

from functions import *

#class handling the normalized data

#dictionary from file name to content

#function to build classifiers for the datasets in dict


class NormalizedDataHandler:
    
    def __init__(self, file_list):
        data_dict = dict()
        for file in file_list:
            data_name = os.path.basename(os.path.splitext(file)[0])
            data_content = pd.read_csv(file, sep = '\t')
            data_dict[data_name] = data_content
        self.data = data_dict 

    #private classification metods running on SINGLE method
    def __dt_classification(self, data_name):
        data = self.data[data_name]
        
        #organize data
        d_pivot=data.pivot(index='sample_id', columns='ab_id', values='ab_count_normalized')
        d_organized = pd.DataFrame(d_pivot.to_records())
        feature_array = d_organized.drop('sample_id',1).values.tolist()

        #assert no nans imported
        array_sum = np.sum(feature_array)
        if(np.isnan(array_sum)):
            print('ERROR: Detected nan value in dataset: ABORT!')
            exit(1)

        colum_order_from_feature_array = d_organized.loc[:,'sample_id']
        sample_cluster_d = data.loc[:,['sample_id','cluster_id']]
        sample_cluster_d = (sample_cluster_d.drop_duplicates(subset=['sample_id']))
        sample_cluster_d['sample_id'] = pd.Categorical(sample_cluster_d['sample_id'], colum_order_from_feature_array)
        sample_cluster_d = sample_cluster_d.sort_values("sample_id")
        class_array = sample_cluster_d.loc[:, 'cluster_id'].values

        X=feature_array
        y=class_array

        #run classification
        cv_inner = KFold(n_splits=10, shuffle=True, random_state=1)
        model = DecisionTreeClassifier(random_state=1)
        params = {"max_depth": [3, None],
              "max_features": range(10, len(feature_array[1])),
              "min_samples_leaf": range(1, 30),
              "criterion": ["gini", "entropy"]}
        # define search
        search = RandomizedSearchCV(model, params, n_iter = 5, scoring='accuracy', n_jobs=1, cv=cv_inner, refit=True)
        # configure the cross-validation procedure
        cv_outer = KFold(n_splits=5, shuffle=True, random_state=1)
        # execute the nested cross-validation
        scores = cross_val_score(search, X, y, scoring='accuracy', cv=cv_outer, n_jobs=-1)
        # report performance
        print('Accuracy[%s] : %.3f (%.3f)' % (data_name, np.mean(scores), np.std(scores)))
        
    #public classification metods running on ALL normalization methods
    #def knn_clasification():

    #def lg_classification():

    def dt_classification(self):
        for key in self.data:
            self.__dt_classification(key)
