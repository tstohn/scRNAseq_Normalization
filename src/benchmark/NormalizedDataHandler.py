from operator import index
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

import shutil
#from cv2 import rotate
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as pltCol
from random import randint
from scipy import stats
import itertools
import math
import random
import collections
from ctypes import *
import sys
sys.path.insert(1, 'src/methods/KnnSimilarity')
from KnnSimilarity import calculate_knn_overlap
from collections import defaultdict

import sklearn
from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score, cross_val_predict, KFold, GridSearchCV, RandomizedSearchCV
from sklearn.tree import *
from sklearn.neighbors import KNeighborsClassifier
from sklearn.manifold import TSNE
from sklearn.metrics import confusion_matrix
from scipy.stats import zscore
import graphviz
import seaborn as sns

import threading

from alibi_detect.cd import MMDDrift #MMD

import sys
sys.path.append('./src/methods/ToolBox')
from functions import *

#class handling the normalized data (classifies, tsne visualisation)
class NormalizedDataHandler:

    def __init_classification(self, data_name, classification_dict):
        data = self.data[data_name]
        
        #dictionary where key is the dataset (Groundtruth, TMM, CLR_COMPOSITIONS,..) and value a dataframe of corelations (see ab_spearman_correlation)
        #dataframe has form index, SPvalues, Pvalues, AB1, AB2, cluster where index is of form AB1_AB2 
        self.correlations = {}

        #organize data
        d_pivot=data.pivot(index='sample_id', columns='ab_id', values='ab_count_normalized')
        d_organized = pd.DataFrame(d_pivot.to_records())
        feature_array = d_organized.drop('sample_id', axis=1).values.tolist()
        feature_ids = d_organized.drop('sample_id', axis=1).columns.tolist()

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
    
    def __init__(self, file_list, groundtruth, deleteBenchmark, threads):
        data_dict = dict()
        classification_dict = dict()
        #threads used for SKLEARN
        self.threads = threads

        #generate dictionary of normalized data
        self.groundtruth = pd.DataFrame()
        for file in file_list:
            data_name = os.path.basename(os.path.splitext(file)[0])
            if("Groundtruth" in data_name):
                data_content = pd.read_csv(file, sep = '\t')
                self.groundtruth = data_content
            data_content = pd.read_csv(file, sep = '\t')
            data_dict[data_name] = data_content

        #dictionary with data for all normalization methods
        sortedDataDict = dict( sorted(data_dict.items(), key=lambda x: x[0].lower()) )
        self.data = sortedDataDict
        self.class_data = classification_dict
        #generate dictionary of classification arrays for normalized data
        for key in self.data:
            self.__init_classification(key, classification_dict)

        #open output file for writing
        folder_path = os.path.dirname(file_list[0])
        folder_name = os.path.basename(folder_path)

        if not os.path.exists("bin/BENCHMARKED_DATASETS/"+folder_name):
            os.mkdir("bin/BENCHMARKED_DATASETS/"+folder_name)
        elif(os.path.exists("bin/BENCHMARKED_DATASETS/"+folder_name) & deleteBenchmark):
            shutil.rmtree("bin/BENCHMARKED_DATASETS/"+folder_name)
            os.mkdir("bin/BENCHMARKED_DATASETS/"+folder_name)

        if not os.path.exists("bin/BENCHMARKED_DATASETS/"+folder_name + "/Overview"):
            os.mkdir("bin/BENCHMARKED_DATASETS/"+folder_name + "/Overview")

        if not os.path.exists("bin/BENCHMARKED_DATASETS/"+folder_name + "/Classification"):
            os.mkdir("bin/BENCHMARKED_DATASETS/"+folder_name + "/Classification")
        if not os.path.exists("bin/BENCHMARKED_DATASETS/"+folder_name + "/Classification/ConfusionMatrices"):
            os.mkdir("bin/BENCHMARKED_DATASETS/"+folder_name + "/Classification/ConfusionMatrices")

        if not os.path.exists("bin/BENCHMARKED_DATASETS/"+folder_name + "/MMDMatrices"):
            os.mkdir("bin/BENCHMARKED_DATASETS/"+folder_name + "/MMDMatrices")

        if not os.path.exists("bin/BENCHMARKED_DATASETS/"+folder_name + "/TSNE_ClusterVisualization"):
            os.mkdir("bin/BENCHMARKED_DATASETS/"+folder_name + "/TSNE_ClusterVisualization")
        if not os.path.exists("bin/BENCHMARKED_DATASETS/"+folder_name + "/SpearmanCorrelations"):
            os.mkdir("bin/BENCHMARKED_DATASETS/"+folder_name + "/SpearmanCorrelations")
        if not os.path.exists("bin/BENCHMARKED_DATASETS/"+folder_name + "/Results"):
            os.mkdir("bin/BENCHMARKED_DATASETS/"+folder_name + "/Results")
        
        #result files
        
        #files origionally used for treatment classification/ now replaced with knn for clusterability
        #self.results = open("bin/BENCHMARKED_DATASETS/"+folder_name+"/Results/treatmentAccuracy.tsv", "w+")
        #self.results.write("NORMALIZATION_METHOD" + "\t" + "CLASSIFICATION_METHOD" + "\t" + "ACCURACY_MEAN" + "\t" + "ACCURACY_SD" + "\n")

        self.sp_results = open("bin/BENCHMARKED_DATASETS/"+folder_name+"/Results/spearmanCorrelations.tsv", "w+")
        self.sp_results.write("NORMALIZATION_METHOD" + "\t" + "SPEARMAN_CORRELATION_MEAN" + "\t" + "SPEARMAN_PVALUE_MEAN" + "\n")

        self.rmsd_sp_results = open("bin/BENCHMARKED_DATASETS/"+folder_name+"/Results/spearmanRMSD.tsv", "w+")
        self.rmsd_sp_results.write("CLUSTER_ID" + "\t" + "NORMALIZATION_METHOD" + "\t" + "SPEARMAN_CORRELATIONS_RMSD" + "\n")

        self.abCorr = open("bin/BENCHMARKED_DATASETS/"+folder_name+"/Results/ABSpearmanCoeff.tsv", "w+")
        self.abCorr.write("CLUSTER_ID" + "\t" + "NORMALIZATION_METHOD" + "\t" + "AB1" + "\t" + "AB2" + "\t" + "GroundtruthCorr" "\t" + "NormCorr\n")

        self.corrDetection = open("bin/BENCHMARKED_DATASETS/"+folder_name+"/Results/PercentageCorrelationDetection.tsv", "w+")
        self.corrDetection.write("cluster"+ "\t" + "norm"+ "\t" + "minCorr"+ "\t" + "TP"+ "\t" + "TN"+ "\t" + "FP"+ "\t" + "FN" + "\n")

        #header is written when stored
        self.knnOverlapFilePath = "bin/BENCHMARKED_DATASETS/"+folder_name+"/Results/knnOverlap.tsv"
        self.clusteringFilePath = "bin/BENCHMARKED_DATASETS/"+folder_name+"/Results/Clustering.tsv"

        self.dataset_name = folder_name
        self.folder_path = ("bin/BENCHMARKED_DATASETS/"+folder_name+"/")

    def __calculate_scores(self, model, params, X, y, global_scores, lock, random_factor=None):
        cv_inner = KFold(n_splits=10, shuffle=True, random_state=random_factor)
        #do seach non threaded, but thread only in cross_val, Linux server has not the dependancies to disentangle those and assign right thread numbers
        search = RandomizedSearchCV(model, params, n_iter = 30, scoring='accuracy', n_jobs=1, cv=cv_inner, refit=True, random_state=random_factor)
        cv_outer = KFold(n_splits=5, shuffle=True, random_state=random_factor)
        scores = cross_val_score(search, X, y, scoring='accuracy', cv=cv_outer, n_jobs=self.threads)
        lock.acquire()
        for s in scores:
            global_scores.append(s)
        lock.release()

    def __classify(self, data_name, model, params, method_string):
        X=self.class_data.get(data_name, {}).get('X')
        y=self.class_data.get(data_name, {}).get('Y')

        threads = []
        global_scores = []
        lock = threading.Lock()
        for t in range(1):
            x = threading.Thread(target=self.__calculate_scores, args=(model, params, X, y, global_scores, lock))
            threads.append(x)
            x.start()
        for x in threads:
            x.join()

        cv_inner = KFold(n_splits=10, shuffle=True, random_state=1)
        search = RandomizedSearchCV(model, params, n_iter = 30, scoring='accuracy', n_jobs=1, cv=cv_inner, refit=True)
        cv_outer = KFold(n_splits=5, shuffle=True, random_state=1)
        scores = cross_val_score(search, X, y, scoring='accuracy', cv=cv_outer, n_jobs=self.threads)
        #plot a confusion matrix
        y_pred = cross_val_predict(search, X, y, cv=cv_outer, n_jobs=self.threads)
        label_list = np.unique(y_pred)
        conf_mat = confusion_matrix(y, y_pred, labels =label_list)
        df_cm = pd.DataFrame(conf_mat, range(len(label_list)), range(len(label_list)))
        plt.figure(figsize=(10,7))
        sns.set(font_scale=1.4) # for label size
        ax= plt.subplot()
        sns.heatmap(df_cm, annot=True, annot_kws={"size": 16}) # font size
        # labels, title and ticks
        ax.set_xlabel('Predicted labels');ax.set_ylabel('True labels')
        ax.set_title('Confusion Matrix')
        ax.xaxis.set_ticklabels(label_list); ax.yaxis.set_ticklabels(label_list)
        plt.savefig(self.folder_path + "/Classification/ConfusionMatrices/" + data_name + method_string + "_CFMatrix.png", dpi=199)
        plt.close()

        printToTerminalOnce('%s => %s Accuracy[%s] : %.3f (%.3f)' % (self.dataset_name, method_string, data_name, np.mean(global_scores), np.std(global_scores)))
        self.results.write(data_name + "\t" + method_string + "\t" + str(round(np.mean(scores), 2)) + "\t" + str(round(np.std(scores), 2)) + "\n")

        return([np.mean(scores),np.std(scores)])

    #private classification metods running on SINGLE method
    #with inner and outer cross validation to get best parameters
    def __knn_optimized_classification(self, data_name, method_string):
        model = KNeighborsClassifier(algorithm='auto')
        feature_array=self.class_data.get(data_name, {}).get('X')
        params = {"n_neighbors": range(1, 30),
                  "leaf_size": range(1,len(feature_array[1])),
                  "p": [1,2],
                  "weights": ["uniform", "distance"],
                  "metric": ["minkowski", "chebyshev"]}
        return(self.__classify(data_name, model, params, method_string))
    
    #knn classification with fixed K and several data-splits
    def __knn_simple_classification(self, key,  knnOverlap, knnMetric, 
                                    test_indices, train_indices,
                                    zscoreScaling = True):

        dataTmp = self.data[key]
        dataTmp = dataTmp[["cluster_id", "sample_id", "ab_id", "ab_count_normalized"]]

        # Pivot the DataFrame to a wider format
        data = dataTmp.pivot(index='sample_id', columns='ab_id', values='ab_count_normalized').reset_index()
        data = data.sort_index()
                
        columns_to_normalize = [col for col in data.columns if col.startswith('AB')]
        if(zscoreScaling):
            data[columns_to_normalize] = data[columns_to_normalize].apply(zscore)

        classLabels = dataTmp[['sample_id', 'cluster_id']].drop_duplicates()
        data = data.merge(classLabels, on='sample_id')

        X = data.drop(['sample_id', 'cluster_id'], axis=1)
        y = data['cluster_id']

        # Split dataset into train and test sets
        X_train = X.iloc[train_indices]
        X_test = X.iloc[test_indices]
        y_train = y.iloc[train_indices]
        y_test = y.iloc[test_indices]
        
        # Create KNN classifier with fixed parameters
        knn = KNeighborsClassifier(n_neighbors=knnOverlap, metric=knnMetric)
        knn.fit(X_train, y_train)
        
        # Make predictions on the test set
        y_pred = knn.predict(X_test)
        # Calculate accuracy
        accuracy = sklearn.metrics.accuracy_score(y_test, y_pred)

        return([accuracy, y_test.tolist(), y_pred.tolist()])

    def __dt_classification(self, data_name, method_string):
        model = DecisionTreeClassifier(random_state=1)
        feature_array=self.class_data.get(data_name, {}).get('X')
        params = {"max_depth": [3, None],
              "max_features": range(10, len(feature_array[1])),
              "min_samples_leaf": range(1, 30),
              "criterion": ["gini", "entropy"]}
        return(self.__classify(data_name, model, params, method_string))

    def calculate_MMD_between_treatments(self):
        #significance test for between treatment difference 
        for key in self.data:
            data = self.data[key]

            #MMD difference between treatments
            #for every two subsets in cluster_id column calculate their MMDDrift
            cluster_ids = data.cluster_id.unique()
            heat_map = pd.DataFrame(0, columns = cluster_ids, index = cluster_ids)
            signi_map = pd.DataFrame(0, columns = cluster_ids, index = cluster_ids)
            for cluster_idx_1 in range(len(cluster_ids)-1):
                data_1 = data[data["cluster_id"]==cluster_ids[cluster_idx_1]]
                mmd_data_prefiltered = data_1.loc[:,['sample_id','ab_id', 'ab_count_normalized']]
                mmd_data_pivotted = mmd_data_prefiltered.pivot(index = "sample_id", columns='ab_id', values='ab_count_normalized')
                mmd_data_1 = mmd_data_pivotted.values
                try:
                    cd = MMDDrift(mmd_data_1, backend='tensorflow', p_val=.05)
                    for cluster_idx_2 in range(cluster_idx_1+1, len(cluster_ids)):
                        data_2 = data[data["cluster_id"]==cluster_ids[cluster_idx_2]]
                        mmd_data_prefiltered_2 = data_2.loc[:,['sample_id','ab_id', 'ab_count_normalized']]
                        mmd_data_pivotted_2 = mmd_data_prefiltered_2.pivot(index = "sample_id", columns='ab_id', values='ab_count_normalized')
                        mmd_data_2 = mmd_data_pivotted_2.values

                        mmd = cd.predict(mmd_data_2, return_p_val=True, return_distance=True)
                        heat_map.loc[cluster_ids[cluster_idx_1], cluster_ids[cluster_idx_2]] = mmd["data"]['distance']
                        signi_map.loc[cluster_ids[cluster_idx_1], cluster_ids[cluster_idx_2]]  =mmd["data"]['is_drift']
                except:
                    continue
            #visualize as heatmap
            plt.figure(figsize=(10,7))
            sns.set(font_scale=1.4) # for label size
            ax= plt.subplot()
            sns.heatmap(heat_map, mask = signi_map == 0, cbar=False, annot=True, annot_kws={"size": 18, "weight": "bold"})
            # labels, title and ticks
            ax.set_xlabel('clusters')
            ax.set_ylabel('clusters')
            ax.set_title('MMD Matrix')
            #ax.xaxis.set_ticklabels(label_list)
            #ax.yaxis.set_ticklabels(label_list)
            plt.savefig(self.folder_path + "MMDMatrices/" + key + ".png", dpi=199)
            plt.close()

    def __draw_tsne(self, data_name):
        X=self.class_data.get(data_name, {}).get('X')
        Y=self.class_data.get(data_name, {}).get('Y')
        tsne = TSNE(n_components=2, verbose=0, perplexity=30, n_iter=400).fit_transform(np.array(X))

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
        plt.savefig(self.folder_path + "TSNE_ClusterVisualization/" + data_name + "_output.png", dpi=199)
        plt.close()

    def draw_treatment_barplot_for_knn_and_dt(self, knnData, dtData):
        # width of the bars
        barWidth = 0.3
        knn = knnData.iloc[0]
        dt = dtData.iloc[0]
        
        # Choose the height of the error bars (bars1)
        knnVar = knnData.iloc[1]
        dtVar = dtData.iloc[1]
        
        # The x position of bars
        r1 = np.arange(len(knn))
        r2 = [x + barWidth for x in r1]
        
        # Create knn bars
        plt.bar(r1, knn, width = barWidth, color = 'blue', edgecolor = 'black', yerr=knnVar, capsize=7, label='KNN')
        # Create dt bars
        plt.bar(r2, dt, width = barWidth, color = 'yellow', edgecolor = 'black', yerr=dtVar, capsize=7, label='DT')
        
        # general layout
        plt.xticks([r + barWidth for r in range(len(knn))], knnData.columns, rotation = 45, ha='right', rotation_mode='anchor', fontsize = 8)
        plt.yticks(fontsize = 8)
        plt.ylabel('Accuracy', fontsize = 9)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 8)
        plt.tight_layout()
        # Show graphic
        plt.savefig(self.folder_path + "Overview/TreatmentClassification.png", dpi=199)
        plt.close()

    #public classification metods running on ALL normalization methods
    def knn_clasification(self,  knnOverlap, knnMetric):
        result = defaultdict(list)
        y_true = defaultdict(list)
        y_pred = defaultdict(list)
        
        print("KNN CLASIIFICATION GOING ON\n")
        #get some random indices for train/ test set. To use exactly the same idx for all normalized datasets for classification
        groundtruthDf = self.data["Groundtruth"]
        unique_samples = groundtruthDf['sample_id'].unique().tolist()
        total_samples = len(unique_samples)
        # Define the number of samples for the test set
        test_size = 0.2
        test_samples = round(test_size * total_samples)

        n_iterations = 50
        # Repeat the process 50 times
        for i in range(n_iterations):
            #create a test/ train split
            random_indices = list(range(total_samples))
            random.shuffle(random_indices)

            test_indices = random_indices[:test_samples]
            # The rest of the indices are for the training set
            train_indices = random_indices[test_samples:]

            for key in self.data:
                classification = self.__knn_simple_classification(key,  knnOverlap, knnMetric, test_indices, train_indices)
                result[key].append(classification[0])
                y_true[key].append(classification[1])
                y_pred[key].append(classification[2])
        
        y_true_list = {key: [item for sublist in value for item in sublist] for key, value in y_true.items()}
        y_pred_list = {key: [item for sublist in value for item in sublist] for key, value in y_pred.items()}

        return([result, y_true_list, y_pred_list])

    def dt_classification(self):
        result = pd.DataFrame()
        for key in self.data:
            result[key] = self.__dt_classification(key, "DecisionTree")
        return(result)
    
    def plot_classification_boxplots(self, results):
        
        #pivot longer: format is one column per normalization method with all accuracies in rows
        resultPivotted = pd.melt(results, var_name='Normalization', value_name='Accuracy')

        # Create boxplot
        plt.figure(figsize=(8, 6))
        sns.boxplot(data=resultPivotted, x='Normalization', y='Accuracy')
        plt.title('Boxplot of Values by Key')
        plt.xlabel('Normalization method')
        plt.ylabel('Accuracy')        
        plt.savefig(self.folder_path +  "/Classification/Clusterability.png", dpi=199)
        plt.close()
        
        
    def sort_classes(self, elem):
        if elem == 'control':
            return -1 
        else:
            return int(elem.split('_')[-1]) 

    def plot_classification_confusion_matrix(self, y_true, y_pred):
        
        #one confusion matrix per normalization method
        for key in self.data:  
            
            #get unique labels, and move control to front
            #only to give confusion matrix an order
            labels = set(y_true[key])
            sorted_labels = sorted(labels, key=self.sort_classes)
 
            cm= confusion_matrix(y_true[key], y_pred[key], sorted_labels)
            
            plt.figure(figsize=(5, 5))
            sns.heatmap(cm, annot = True, linewidths= 0.5, linecolor="red", fmt=".0f")

            plt.xlabel('Predicted labels')
            plt.ylabel('True labels')
            plt.title('Confusion Matrix')
            # Set tick labels
            plt.xticks(ticks=range(len(sorted_labels)), labels=sorted_labels)
            plt.yticks(ticks=range(len(sorted_labels)), labels=sorted_labels)
            # Save the plot
            plt.savefig(self.folder_path +  "/Classification/ConfusionMatrices/ConfusionKNN" + key + ".png", dpi=199)
            plt.close()
        
    def validate_clusterability(self, knnOverlap, knnMetric):
        #we return the accuracies per norm method, as well as predicted, true labels for conf matrices
        knnResults = self.knn_clasification(knnOverlap, knnMetric)

        #accuracies boxplot
        knnScores = knnResults[0]
        knnScores = pd.DataFrame(knnScores)
        
        knnScores.sort_index(axis=1, inplace=True)
        self.plot_classification_boxplots(knnScores)
        
        knnScores["KnnThreshold"] = knnOverlap
        knnScores["knnMetric"] = knnMetric
        knnScores.to_csv(self.clusteringFilePath, sep='\t', index=False)

        #confusion matrix
        self.plot_classification_confusion_matrix(knnResults[1], knnResults[2])

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
                X_train, X_test = np.array(X)[train_index], np.array(X)[test_index]
                y_train, Y_test = np.array(Y)[train_index], np.array(Y)[test_index]
            #model = model.fit(X_train, y_train)
            cv_inner = KFold(n_splits=5, shuffle=True, random_state=1)
            search = RandomizedSearchCV(model, params, n_iter = 20, scoring='accuracy', n_jobs=self.threads, cv=cv_inner, refit=True)
            search = search.fit(X_train, y_train)
            model = search.best_estimator_
            dot_data = export_graphviz(model, out_file=None, 
                                    feature_names=self.class_data.get(data_name, {}).get('FEATURES'),  
                                    class_names=np.unique(y_train),
                                    filled=True,
                                    precision = 4)
            graph = graphviz.Source(dot_data, format="png") 
            graph.render(data_name + "decision_tree_graphivz")

    def __calculate_all_spearman_correlations(self, data):
        
        data_subset = data.loc[:,['sample_id','ab_id', 'ab_count_normalized']]
        d_pivot=data_subset.pivot(index = "sample_id", columns='ab_id', values='ab_count_normalized')
        columns = d_pivot.columns.tolist()

        #dataframe to fill with spearkman values
        column_names = ["index", "SPvalues", "Pvalues", "AB1", "AB2"]
        correlations = pd.DataFrame(columns = column_names)
        for col_a, col_b in itertools.combinations(columns, 2):

            SPvalues, Pvalues = stats.spearmanr(d_pivot.loc[:, col_a], d_pivot.loc[:, col_b])
            newIndex = col_a + '_' + col_b

            newCorrLine = pd.DataFrame.from_dict({'index': [newIndex], 'SPvalues': [SPvalues], 'Pvalues': [Pvalues], 'AB1': [col_a], 'AB2': [col_b]}, orient='columns')
            correlations = pd.concat([correlations, newCorrLine], ignore_index = True)
            
        return(correlations)
    
    def __calculate_all_spearman_correlations_per_cluster(self, data):
        
        dataGroundtruth = self.data["Groundtruth"].copy(deep = True)
        correlationsPerCluster = {}
        for cluster in dataGroundtruth['cluster_id'].unique():
            data = data[data["cluster_id"] == cluster]

            data_subset = data.loc[:,['sample_id','ab_id', 'ab_count_normalized']]
            d_pivot=data_subset.pivot(index = "sample_id", columns='ab_id', values='ab_count_normalized')
            columns = d_pivot.columns.tolist()

            #dataframe to fill with spearkman values
            column_names = ["index", "SPvalues", "Pvalues", "AB1", "AB2"]
            correlations = pd.DataFrame(columns = column_names)
            for col_a, col_b in itertools.combinations(columns, 2):

                SPvalues, Pvalues = stats.spearmanr(d_pivot.loc[:, col_a], d_pivot.loc[:, col_b])
                newIndex = col_a + '_' + col_b

                newCorrLine = pd.DataFrame.from_dict({'index': [newIndex], 'SPvalues': [SPvalues], 'Pvalues': [Pvalues], 'AB1': [col_a], 'AB2': [col_b]}, orient='columns')
                correlations = pd.concat([correlations, newCorrLine], ignore_index = True)
                
            correlationsPerCluster[cluster] = correlations
        return(correlationsPerCluster)

    def save_all_correlations(self, data, key, cluster):
        for row in data.itertuples():
            self.abCorr.write(cluster + "\t" + key + "\t" + row.AB1 + "\t" + row.AB2 + "\t" + str(round(row.SPvalues, 4)) + "\t" + str(round(row.SPvalues_norm, 4)) + "\n")

    def ab_spearman_correlation_graph(self, filter = ""):
        plt.figure(figsize=(16,10))
        lengend_labels = []

        for key in self.data:
            data = self.data[key].copy(deep = True)

            #speacial line for a data set where phospho proteins were excluded due to expected correlations
            if(filter != ""):
                data = data[data['ab_type'] == filter]
            spearmanValues = self.__calculate_all_spearman_correlations(data.copy(deep = True))
            sp_mean = np.mean(spearmanValues.SPvalues)
            p_mean = np.mean(spearmanValues.Pvalues)

            self.sp_results.write(key + "\t"+ str(round(sp_mean,4)) + "\t" + str(round(p_mean,4)) + "\n")
            sns.distplot(x=spearmanValues.SPvalues, hist=False, kde=True)
            lengend_labels.append(key)

        plt.legend(labels=lengend_labels)
        plt.savefig(self.folder_path +  "SpearmanCorrelations/spearman_correlations_global_" + filter + ".png", dpi=199)
        plt.close()
        
    def calculate_all_correlations(self):
        
        for key in self.data:
            dataAllClusters = self.data[key].copy()
            column_names = ["cluster","index", "SPvalues", "Pvalues", "AB1", "AB2"]
            correlations = pd.DataFrame(columns = column_names)
            for cluster in dataAllClusters['cluster_id'].unique():
                data = dataAllClusters[dataAllClusters["cluster_id"] == cluster]
                correlationsTmp = self.__calculate_all_spearman_correlations(data.copy(deep = True))
                correlationsTmp["cluster"] = cluster
                correlations = pd.concat([correlations, correlationsTmp], ignore_index=True)

            self.correlations[key] = correlations
            
    def ab_correlation_evaluation(self):

        #build a dataFrame that calculates all correlations between proteins for grondTruth
        groundTruthCorrelations = self.correlations["Groundtruth"]
        
        normRMSDData = {}
        #for every norm method
        for key in self.data:
            if("Groundtruth" in key):
                continue
            #get cluster specific correlations
            dataAllClusters = self.data[key].copy(deep = True)
            for cluster in dataAllClusters['cluster_id'].unique():

                #calculate same dataFrame and add GroundTruth to it
                normCorrelations = self.correlations[key]
                normCorrelations = normCorrelations[normCorrelations["cluster"] == cluster]
                normCorrelations.columns = ["cluster", "index","SPvalues_norm", "Pvalues_norm", "AB1_norm", "AB2_norm"]
                
                groundTruthCorrelationsTmp = groundTruthCorrelations[groundTruthCorrelations["cluster"] == cluster]
                
                #dropping nan values, important for e.g. subsampling
                groundTruthCorrelationsTmp.dropna(subset = ["SPvalues"], inplace=True)
                normCorrelations.dropna(subset = ["SPvalues_norm"], inplace=True)

                result = pd.merge(groundTruthCorrelationsTmp, normCorrelations, how="inner", on = "index")

                #calculate the RMSD and safe it to RMSD matrix
                mse = sklearn.metrics.mean_squared_error(result.SPvalues, result.SPvalues_norm)
                rmsd = math.sqrt(mse)
                normRMSDData[key] = rmsd
                self.rmsd_sp_results.write(cluster + "\t" + key + "\t"+ str(round(rmsd,4)) + "\n")

                #wite all correlations to file for AB pairs
                self.save_all_correlations(result, key, cluster)

        #draw barplot of all RMSDs
        values = normRMSDData.values()
        normMethods = normRMSDData.keys()
        y_pos = np.arange(len(normMethods))
        plt.bar(y_pos, values, color = 'blue', edgecolor = 'black')
        plt.ylabel('RMSD of correlations between AB counts', fontsize = 10)
        plt.xticks(y_pos, normMethods, rotation = 45, ha='right', rotation_mode='anchor', fontsize = 8)
        plt.yticks(fontsize = 8)
        plt.tight_layout()
        plt.savefig(self.folder_path + "Overview/CorrelationRMSD.png", dpi=199)
        plt.close()

    def ab_spearman_correlation(self, groundtruth = False, filter = ""):
        #plot graph of spearman correlations
        self.ab_spearman_correlation_graph(filter)

        #generate barplot of RMSD for normalized correlations and groundtruth correlations
        if( (not filter) & (groundtruth) ):
            self.ab_correlation_evaluation()

    def __calc_fold_changes_between_totalABCounts_and_ABmin(self, data):

        summedData = data.groupby(['sample_id'])['ab_count_normalized'].sum().reset_index()
        summedData.rename(columns={'ab_count_normalized' :'ab_count_total'}, inplace=True)

        summedData = pd.merge(summedData, data, how="outer", on = ["sample_id"])
        summedData.drop_duplicates(inplace=True)
        
        data_subset = summedData.loc[:,['sample_id', 'ab_count_total']]
        totalCountDict = dict(zip(data_subset.sample_id, data_subset.ab_count_total))     

        minCount = min(totalCountDict.values())   
        foldChanges = {}
        for sample in totalCountDict.keys():
            foldChanges[sample] = totalCountDict[sample]/minCount

        result = pd.DataFrame.from_dict(foldChanges, orient='index')
        result.columns = ['foldChange']
        result.reset_index(inplace=True)

        return(result)

    def validate_normalizedData_against_groundTruth(self):
        
        normRMSDData = {}
        dataTruth = self.groundtruth.copy(deep = True)
        dataTruth = dataTruth[['sample_id', 'ab_id', 'ab_count_normalized']]
        dataTruth = self.__calc_fold_changes_between_totalABCounts_and_ABmin(dataTruth)
        dataTruth.rename(columns={'foldChange' :'foldChange_truth'},inplace=True)
        #for every norm method
        for key in self.data:
            if(key == "CLR"):
                continue
            if("Groundtruth" in key):
                continue
            dataNorm = self.data[key].copy(deep = True)
            dataNorm = dataNorm[['sample_id', 'ab_id', 'ab_count_normalized']]
            dataNorm = self.__calc_fold_changes_between_totalABCounts_and_ABmin(dataNorm)
            dataNorm.rename(columns={'foldChange' :'foldChange_norm'},inplace=True)

            result = pd.merge(dataNorm, dataTruth, how="outer", on = ["index"])

            #calculate the RMSD and safe it to RMSD matrix
            mse = sklearn.metrics.mean_squared_error(result.foldChange_norm, result.foldChange_truth)
            rmsd = math.sqrt(mse)
            normRMSDData[key] = rmsd

        #draw barplot of all RMSDs
        values = normRMSDData.values()
        normMethods = normRMSDData.keys()
        y_pos = np.arange(len(normMethods))
        plt.bar(y_pos, values, color = 'blue', edgecolor = 'black')
        plt.ylabel('RMSD of the ration \n(total AB count/min(total AB count)) per sample', fontsize = 10)
        plt.xticks(y_pos, normMethods, rotation = 45, ha='right', rotation_mode='anchor', fontsize = 8)
        plt.yticks(fontsize = 8)
        plt.tight_layout()
        plt.savefig(self.folder_path + "Overview/ABCountRMSD.png", dpi=199)
        plt.close()

    def validate_correlations(self, corrCutoff = 0.5, stepSize = 0.05):
        
        #get list of all indeces (AB1_AB2) of TP correlations
        column_names = ["cluster", "norm", "minCorr", "TP", "TN", "FP", "FN"]
        correlationRocResults = pd.DataFrame(columns = column_names)
        groundTruthCorrelations = self.correlations["Groundtruth"].copy(deep = True)
        
        groundtruthData = self.data["Groundtruth"].copy(deep = True)
        clusters = groundtruthData['cluster_id'].unique()
        truePosCorrelations = {key: [] for key in clusters}

        for cluster in clusters:
            groundTruthTmp = groundTruthCorrelations[groundTruthCorrelations["cluster"] == cluster]
            for index, row in groundTruthTmp.iterrows():
                if( abs(row["SPvalues"]) > corrCutoff):
                    truePosCorrelations[cluster].append(row["index"])
                    
        #for all norm methods        
        for key in self.data:
            if("Groundtruth" in key):
                continue
            
            #for all clusters
            for cluster in clusters:
                normData = self.correlations[key].copy(deep = True)
                normData = normData[normData["cluster"] == cluster]
                
                cutOffTmp = 1.0
                while 0.0 <= cutOffTmp:
                    
                    #get correlastions > threshold
                    detectedCorrelations = []
                    for index, row in normData.iterrows():
                        if( abs(row["SPvalues"]) > cutOffTmp):
                            detectedCorrelations.append(row["index"])
                           
                    #count number of TP,TN,FP,FN indices and add to result
                    tp = sum(corr in truePosCorrelations[cluster] for corr in detectedCorrelations)
                    fp = len(detectedCorrelations) - tp
                    
                    tn = len(normData.index) - len(truePosCorrelations[cluster]) - fp
                    fn = len(normData.index) - tn - len(detectedCorrelations)
                            
                    data = {"cluster": [cluster],"norm": [key], "minCorr": [cutOffTmp],
                            "TP":[tp], "TN":[tn], "FP":[fp], "FN":[fn]}
                    resultTmp = pd.DataFrame(data)
                    correlationRocResults = pd.concat([correlationRocResults, resultTmp], ignore_index=True)

                    cutOffTmp = cutOffTmp - stepSize
                    self.corrDetection.write(cluster + "\t" + key + "\t"+ str(round(cutOffTmp, 2))+ "\t" + str(tp) + "\t" + str(tn) + "\t" + str(fp) + "\t" + str(fn) + "\n")

    #mean values are in the end divided by the smallest AB value
    def __get_AB_mean_dataFrame(self, data):

        batchIds = data["batch_id"].unique()
        result = pd.DataFrame()
        for batch in batchIds:
            dataBatch = data[data.batch_id == batch]
            dataBatchMeanAbCounts = dataBatch.groupby(['ab_id'])['ab_count_normalized'].mean().reset_index()
            dataBatchMeanAbCounts.rename(columns={'ab_count_normalized' :'ab_count_mean'}, inplace=True)
            dataBatchMeanAbCounts = dataBatchMeanAbCounts[['ab_id', 'ab_count_mean']]
            minValue = min(dataBatchMeanAbCounts["ab_count_mean"])
            dataBatchMeanAbCounts["ab_count_mean"] = dataBatchMeanAbCounts["ab_count_mean"]/minValue
            dataBatchMeanAbCounts.set_index('ab_id', inplace = True)
            result[batch] = dataBatchMeanAbCounts['ab_count_mean']

        return(result)

    def validate_batch_effect(self):

        #for every batch
        normMethodsBatchDiff = pd.DataFrame()
        for key in self.data:
            if("Groundtruth" in key):
                continue
            if("Simulation" in key):
                continue
            batchDifferenceDict = {}
            data = self.data[key].copy(deep = True)
            meanData = self.__get_AB_mean_dataFrame(data)
            batches = meanData.columns.to_list()
            for batch_1, batch_2 in itertools.combinations(batches, 2):
                mse = sklearn.metrics.mean_squared_error(meanData[batch_1], meanData[batch_2])
                rmsd = math.sqrt(mse)
                batchDifferenceDict[str(batch_1) + "_" + str(batch_2)] = rmsd            
            result = pd.DataFrame(batchDifferenceDict.items(), columns=['batch_correlation', 'rmsd_value'])
            normMethodsBatchDiff[key] = result["rmsd_value"]

        ax= plt.subplot()
        ax = normMethodsBatchDiff.plot.bar(rot=0)
        # labels, title and ticks
        ax.set_ylabel('RMSD of mean AB count between batches', fontsize = 10)
        ax.set_title('Batch effect', fontsize = 10)
        plt.xticks(rotation = 45, ha='right', rotation_mode='anchor', fontsize = 8)
        plt.yticks(rotation = 45, ha='right', rotation_mode='anchor', fontsize = 6)
        plt.xlabel('RMSD between all combinations of batches', fontsize = 10)

        plt.tight_layout()
        plt.savefig(self.folder_path + "Overview/BatchEffect.png", dpi=199)
        plt.close()
        
    #treatment conditions are enumerated from 0 to X
    def __parseTreatmentConditionInfluence(self, diffExProteins, diffExFactors):
        dict = {}
        conditionId = 0
        for cond in diffExProteins:
            abId = 0
            for abIdStr in cond:
                #conditions are enumerated from 1
                dict[(abIdStr, "cluster_" + str(conditionId))] = diffExFactors[conditionId][abId]
                abId = abId + 1
            conditionId = conditionId + 1

        return(dict)

    def validate_log2foldchange(self, diffExProteins, diffExFactors):

        abIdConditionToFactorMap = self.__parseTreatmentConditionInfluence(diffExProteins, diffExFactors)
        treatmentProts = []
        #get all diff proteins
        for vec in diffExProteins:
            treatmentProts.extend(vec)
        
        fig, axs = plt.subplots(len(self.data), len(treatmentProts),figsize=(17,12))
        i = 0
        for key in self.data:
            data = self.data[key].copy(deep = True)
            data = data[data.ab_id.isin(treatmentProts)]
            j = 0
            for prot in treatmentProts:
                dict = {}
                dataTmp = data[data["ab_id"] == prot]
                for condition in dataTmp.cluster_id.unique():
                    condition = str(condition)
                    dict[condition] = dataTmp.loc[dataTmp["cluster_id"] == condition, "ab_count_normalized"]
                if(not dict):
                    printToTerminalOnce("Error: There is no data for conditions. E.g. AB names do not match anymore (AB duplicates)")
                od = collections.OrderedDict(sorted(dict.items()))
                box_plots = axs[i, j].boxplot(od.values(), patch_artist=True, showfliers=False)
                axs[i, j].set_xticklabels(od.keys(), size = '10')
                axs[i, j].set_title(key + " " + prot, size = '12')
                pos = "green"
                neg = "red"
                non = ""
                colList = []
                for cond in od.keys():
                    if((prot,cond) in abIdConditionToFactorMap):
                        if(abIdConditionToFactorMap[(prot,cond)]>1):
                            colList.append(pos)
                        else:
                            colList.append(neg)
                    else:
                        colList.append(non)
                for b,box in enumerate(box_plots['boxes']):
                    if(colList[b] == ""):
                        color_with_alpha = pltCol.colorConverter.to_rgba("white", 0.0)
                        box.set_facecolor(color_with_alpha)
                        continue
                    box.set_facecolor(colList[b])
                    box.set_alpha(0.7)
                    box.set_edgecolor("black")
                #plt.setp(box1["boxes"], facecolor="red")
                #add jitter
                keyPos = 0
                for key, value in od.items():
                    y = value
                    x = np.random.normal(keyPos + 1, 0.1, len(y))
                    #axs[i, j].scatter(x, y, alpha = 0.4, s = 0.2)
                    keyPos = keyPos + 1
                j = j + 1
            i = i+1
        plt.tight_layout()
        plt.savefig(self.folder_path + "Overview/TreatmentWOOutliers.png", dpi=199)
        plt.close()

    def validate_knn_overlap(self, knnOverlap, knnMetric):
        
        #Steps
        #bring data frames into right format
        #calculate knn overlap between them
        #write overlap into a result table

        #result: row is single cell, column all norm methods, value percent overlap
        normMethods = self.data.keys()
        result = None

        #for every norm method
        for key in self.data:

            if("Groundtruth" in key):
                continue
            
            dataTmp = self.data[key].copy(deep = True)
            groundtruth = self.groundtruth.copy(deep = True)
            
            #if we run a gradient of various KNN thresholds (knnOverlap == 0)
            print(knnOverlap)
            if(knnOverlap == 0):
                #in this case we take neighborhood numbers from 5% of samples to 95%
                sampleNumber = len(dataTmp.sample_id.unique())
                for percentage in np.arange(0.05, 1.0, 0.05):
                    knn = round(percentage * sampleNumber)
                    overlapDist = calculate_knn_overlap(dataTmp, groundtruth, knn, knnMetric)
                    overlapDist['NORM_METHOD'] = str(key)

                    if(result is None):
                        result = overlapDist
                    else:
                        result = pd.concat([result, overlapDist], ignore_index=True)
        
            else: #if we only use a single threshold > 0
                knn = knnOverlap
                overlapDist = calculate_knn_overlap(dataTmp, groundtruth, knn, knnMetric)
                overlapDist['NORM_METHOD'] = str(key)

                if(result is None):
                    result = overlapDist
                else:
                    result = pd.concat([result, overlapDist], ignore_index=True)
        
        #write final results
        result.to_csv(self.knnOverlapFilePath, sep='\t', index=False)





            


