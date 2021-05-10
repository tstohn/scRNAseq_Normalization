import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from random import randint
from scipy import stats

from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score, KFold, GridSearchCV, RandomizedSearchCV
from sklearn.tree import *
from sklearn.neighbors import KNeighborsClassifier
from sklearn.manifold import TSNE

import graphviz
import seaborn as sns

#class handling the normalized data (classifies, tsne visualisation)
class NormalizedDataHandler:

    def __init_classification(self, data_name, classification_dict):
        data = self.data[data_name]

        #organize data
        d_pivot=data.pivot(index='sample_id', columns='ab_id', values='ab_count_normalized')
        d_organized = pd.DataFrame(d_pivot.to_records())
        feature_array = d_organized.drop('sample_id',1).values.tolist()
        feature_ids = d_organized.drop('sample_id',1).columns.tolist()

        #assert no nans imported
        array_sum = np.sum(feature_array)
        assert(not np.isnan(array_sum))

        colum_order_from_feature_array = d_organized.loc[:,'sample_id']
        sample_cluster_d = data.loc[:,['sample_id','cluster_id']]
        sample_cluster_d = (sample_cluster_d.drop_duplicates(subset=['sample_id']))
        sample_cluster_d['sample_id'] = pd.Categorical(sample_cluster_d['sample_id'], colum_order_from_feature_array)
        sample_cluster_d = sample_cluster_d.sort_values("sample_id")
        class_array = sample_cluster_d.loc[:, 'cluster_id'].values

        X=feature_array
        Y=class_array

        classification_dict[data_name] = {"X": X, "Y": Y, "FEATURES": feature_ids}
    
    def __init__(self, file_list):
        data_dict = dict()
        classification_dict = dict()
        #generate dictionary of normalized data
        for file in file_list:
            data_name = os.path.basename(os.path.splitext(file)[0])
            data_content = pd.read_csv(file, sep = '\t')
            data_dict[data_name] = data_content
        #dictionary with data for all normalization methods
        self.data = data_dict
        self.class_data = classification_dict
        #generate dictionary of classification arrays for normalized data
        for key in self.data:
            self.__init_classification(key, classification_dict)

        #open output file for writing
        folder_path = os.path.dirname(file_list[0])
        folder_name = os.path.basename(folder_path)
        if not os.path.exists("bin/BENCHMARKED_DATASETS/"+folder_name):
            os.mkdir("bin/BENCHMARKED_DATASETS/"+folder_name)
        self.results = open("bin/BENCHMARKED_DATASETS/"+folder_name+"/results.tsv", "w+")
        self.results.write("NORMALIZATION_METHOD" + "\t" + "CLASSIFICATION_METHOD" + "\t" + "ACCURACY_MEAN" + "\t" + "ACCURACY_SD" + "\n")

        self.sp_results = open("bin/BENCHMARKED_DATASETS/"+folder_name+"/results_spearman.tsv", "w+")
        self.sp_results.write("NORMALIZATION_METHOD" + "\t" + "SPEARMAN_CORRELATION_MEAN" + "\t" + "SPEARMAN_PVALUE_MEAN" + "\n")

        self.dataset_name = folder_name
        self.folder_path = ("bin/BENCHMARKED_DATASETS/"+folder_name+"/")

    #def __del__(self):
        #self.results.close()

    def __classify(self, data_name, model, params, method_string):
        X=self.class_data.get(data_name, {}).get('X')
        Y=self.class_data.get(data_name, {}).get('Y')

        cv_inner = KFold(n_splits=10, shuffle=True, random_state=1)
        search = RandomizedSearchCV(model, params, n_iter = 30, scoring='accuracy', n_jobs=1, cv=cv_inner, refit=True)
        cv_outer = KFold(n_splits=5, shuffle=True, random_state=1)
        scores = cross_val_score(search, X, Y, scoring='accuracy', cv=cv_outer, n_jobs=-1)
        print('%s => %s Accuracy[%s] : %.3f (%.3f)' % (self.dataset_name, method_string, data_name, np.mean(scores), np.std(scores)))
        
        self.results.write(data_name + "\t" + method_string + "\t" + str(round(np.mean(scores), 2)) + "\t" + str(round(np.std(scores), 2)) + "\n")

    #private classification metods running on SINGLE method
    def __knn_classification(self, data_name, method_string):
        model = KNeighborsClassifier(algorithm='auto')
        feature_array=self.class_data.get(data_name, {}).get('X')
        params = {"n_neighbors": range(1, 30),
                  "leaf_size": range(1,len(feature_array[1])),
                  "p": [1,2],
                  "weights": ["uniform", "distance"],
                  "metric": ["minkowski", "chebyshev"]}
        self.__classify(data_name, model, params, method_string)

    def __dt_classification(self, data_name, method_string):
        model = DecisionTreeClassifier(random_state=1)
        feature_array=self.class_data.get(data_name, {}).get('X')
        params = {"max_depth": [3, None],
              "max_features": range(10, len(feature_array[1])),
              "min_samples_leaf": range(1, 30),
              "criterion": ["gini", "entropy"]}
        self.__classify(data_name, model, params, method_string)

    def __draw_tsne(self, data_name):
        X=self.class_data.get(data_name, {}).get('X')
        Y=self.class_data.get(data_name, {}).get('Y')

        tsne = TSNE(n_components=2, verbose=0, perplexity=30, n_iter=400).fit_transform(X)

        tsne_df = pd.DataFrame({'X':tsne[:,0],
                        'Y':tsne[:,1],
                        'class':Y})

        plt.figure(figsize=(16,10))
        sns.scatterplot(
            x="X", y="Y",
            hue="class",
            #palette=sns.color_palette("viridis", 4),
            data=tsne_df,
            legend="full",
            alpha=1
        )
        plt.savefig(self.folder_path + data_name + "_output.png")

        
    #public classification metods running on ALL normalization methods
    def knn_clasification(self):
        for key in self.data:
            self.__knn_classification(key, "KNN")

    #def lg_classification():

    def dt_classification(self):
        for key in self.data:
            self.__dt_classification(key, "DecisionTree")

    def draw_tsne(self):
        for key in self.data:
            self.__draw_tsne(key)

    def draw_dt_graph(self):
        for data_name in self.data:
            model = DecisionTreeClassifier(random_state=1)
            feature_array=self.class_data.get(data_name, {}).get('X')
            params = {"max_depth": [3, None],
                "max_features": range(10, len(feature_array[1])),
                "min_samples_leaf": range(1, 30),
                "criterion": ["gini", "entropy"]}

            #graph of DC tree
            cv = KFold(n_splits=5, shuffle=True, random_state=1)
            X=self.class_data.get(data_name, {}).get('X')
            Y=self.class_data.get(data_name, {}).get('Y')
            for train_index, test_index in cv.split(X):
                print(train_index)
                X_train, X_test = np.array(X)[train_index], np.array(X)[test_index]
                y_train, Y_test = np.array(Y)[train_index], np.array(Y)[test_index]
            #model = model.fit(X_train, y_train)
            cv_inner = KFold(n_splits=5, shuffle=True, random_state=1)
            search = RandomizedSearchCV(model, params, n_iter = 20, scoring='accuracy', n_jobs=1, cv=cv_inner, refit=True)
            search = search.fit(X_train, y_train)
            model = search.best_estimator_
            dot_data = export_graphviz(model, out_file=None, 
                                    feature_names=self.class_data.get(data_name, {}).get('FEATURES'),  
                                    class_names=np.unique(y_train),
                                    filled=True,
                                    precision = 4)
            graph = graphviz.Source(dot_data, format="png") 
            graph.render(data_name + "decision_tree_graphivz")

    def ab_spearman_correlation(self):
        plt.figure(figsize=(16,10))
        lengend_labels = []
        for key in self.data:
            data = self.data[key]

            #speacial line for a data set where phospho proteins were excluded due to expected correlations
            #data = data[data['ab_type'] == "total"]

            data_subset = data.loc[:,['sample_id','ab_id', 'ab_count_normalized']]
            d_pivot=data_subset.pivot(index = "ab_id", columns='sample_id', values='ab_count_normalized')
            sp, p = stats.spearmanr(d_pivot, axis = 1)
            dim=len(sp[0])
            
            sp_values=list(sp[np.triu_indices(dim,1)])
            p_values=list(p[np.triu_indices(dim,1)])

            sp_mean = np.mean(sp_values)
            p_mean = np.mean(p_values)

            self.sp_results.write(key + "\t"+ str(round(sp_mean,4)) + "\t" + str(round(p_mean,4)) + "\n")
            sns.distplot(x=sp_values, hist=False, kde=True)
            lengend_labels.append(key)
        plt.legend(labels=lengend_labels)
        plt.savefig(self.folder_path +  "spearman_correlations.png")






