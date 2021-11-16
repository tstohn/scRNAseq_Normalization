from operator import index
import os
import shutil
from cv2 import rotate
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

import sklearn
from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score, cross_val_predict, KFold, GridSearchCV, RandomizedSearchCV
from sklearn.tree import *
from sklearn.neighbors import KNeighborsClassifier
from sklearn.manifold import TSNE
from sklearn.metrics import confusion_matrix

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
        groundtruth = None

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
        self.results = open("bin/BENCHMARKED_DATASETS/"+folder_name+"/Classification/results.tsv", "w+")
        self.results.write("NORMALIZATION_METHOD" + "\t" + "CLASSIFICATION_METHOD" + "\t" + "ACCURACY_MEAN" + "\t" + "ACCURACY_SD" + "\n")

        self.sp_results = open("bin/BENCHMARKED_DATASETS/"+folder_name+"/SpearmanCorrelations/results_spearman.tsv", "w+")
        self.sp_results.write("NORMALIZATION_METHOD" + "\t" + "SPEARMAN_CORRELATION_MEAN" + "\t" + "SPEARMAN_PVALUE_MEAN" + "\n")

        self.dataset_name = folder_name
        self.folder_path = ("bin/BENCHMARKED_DATASETS/"+folder_name+"/")

    def __calculate_scores(self, model, params, X, y, global_scores, lock, random_factor=None):
        cv_inner = KFold(n_splits=10, shuffle=True, random_state=random_factor)
        search = RandomizedSearchCV(model, params, n_iter = 30, scoring='accuracy', n_jobs=self.threads, cv=cv_inner, refit=True, random_state=random_factor)
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
        search = RandomizedSearchCV(model, params, n_iter = 30, scoring='accuracy', n_jobs=self.threads, cv=cv_inner, refit=True)
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
    def __knn_classification(self, data_name, method_string):
        model = KNeighborsClassifier(algorithm='auto')
        feature_array=self.class_data.get(data_name, {}).get('X')
        params = {"n_neighbors": range(1, 30),
                  "leaf_size": range(1,len(feature_array[1])),
                  "p": [1,2],
                  "weights": ["uniform", "distance"],
                  "metric": ["minkowski", "chebyshev"]}
        return(self.__classify(data_name, model, params, method_string))

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
    def knn_clasification(self):
        result = pd.DataFrame()
        for key in self.data:
            result[key] = self.__knn_classification(key, "KNN")
        return(result)

    def dt_classification(self):
        result = pd.DataFrame()
        for key in self.data:
            result[key] = self.__dt_classification(key, "DecisionTree")
        return(result)

    def run_treatment_classification(self):
        dtScores = self.dt_classification()
        knnScores = self.knn_clasification()

        dtScores.sort_index(axis=1, inplace=True)
        knnScores.sort_index(axis=1, inplace=True)

        self.draw_treatment_barplot_for_knn_and_dt(knnScores, dtScores)

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
            correlations = {}

            for col_a, col_b in itertools.combinations(columns, 2):
                correlations[col_a + '_' + col_b] = stats.spearmanr(d_pivot.loc[:, col_a], d_pivot.loc[:, col_b])

            result = pd.DataFrame.from_dict(correlations, orient='index')
            result.columns = ['SPvalues', 'Pvalues']
            result.reset_index(inplace=True)

            return(result)

    def ab_spearman_correlation_graph(self, filter = ""):
        plt.figure(figsize=(16,10))
        lengend_labels = []

        for key in self.data:
            print(key)
            data = self.data[key].copy()

            #speacial line for a data set where phospho proteins were excluded due to expected correlations
            if(filter != ""):
                data = data[data['ab_type'] == filter]
            spearmanValues = self.__calculate_all_spearman_correlations(data)
            sp_mean = np.mean(spearmanValues.SPvalues)
            p_mean = np.mean(spearmanValues.Pvalues)

            self.sp_results.write(key + "\t"+ str(round(sp_mean,4)) + "\t" + str(round(p_mean,4)) + "\n")
            sns.distplot(x=spearmanValues.SPvalues, hist=False, kde=True)
            lengend_labels.append(key)

        plt.legend(labels=lengend_labels)
        plt.savefig(self.folder_path +  "SpearmanCorrelations/spearman_correlations_" + filter + ".png", dpi=199)
        plt.close()

    def ab_correlation_evaluation(self):

        #build a dataFrame that calculates all correlations between proteins for grondTruth
        groundTruthCorrelations = self.__calculate_all_spearman_correlations(self.groundtruth.copy())
        
        normRMSDData = {}
        #for every norm method
        for key in self.data:
            if("Groundtruth" in key):
                continue

            data = self.data[key].copy()
            #calculate same dataFrame and add GroundTruth to it
            normCorrelations = self.__calculate_all_spearman_correlations(data)
            normCorrelations.columns = ['index','SPvalues_norm', 'Pvalues_norm']

            #dropping nan values, important for e.g. subsampling
            groundTruthCorrelations.dropna(subset = ["SPvalues"], inplace=True)
            normCorrelations.dropna(subset = ["SPvalues_norm"], inplace=True)

            result = pd.merge(groundTruthCorrelations, normCorrelations, how="inner", on = "index")

            #calculate the RMSD and safe it to RMSD matrix
            mse = sklearn.metrics.mean_squared_error(result.SPvalues, result.SPvalues_norm)
            rmsd = math.sqrt(mse)
            normRMSDData[key] = rmsd
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
        dataTruth = self.groundtruth.copy()
        dataTruth = dataTruth[['sample_id', 'ab_id', 'ab_count_normalized']]
        dataTruth = self.__calc_fold_changes_between_totalABCounts_and_ABmin(dataTruth)
        dataTruth.rename(columns={'foldChange' :'foldChange_truth'},inplace=True)
        #for every norm method
        for key in self.data:
            if(key == "CLR"):
                continue
            if("Groundtruth" in key):
                continue
            dataNorm = self.data[key].copy()
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

    def validate_correlations(self, proteinCorrelations):
        
        #for the proteinCorrelations calcualte the spearman correlation for groundtruth and norm
        indices = []
        columnNames = ["Groundtruth"]
        wantedVarProteinList = []
        for corr in proteinCorrelations:
            prot1 = corr.prot1
            prot2 = corr.prot2
            wantedVarProteinList.append(prot1)
            wantedVarProteinList.append(prot2)
            indices.append(prot1 + '_' + prot2)
        for key in self.data:
            if("Groundtruth" in key):
                continue
            columnNames.append(key)
        correlations = pd.DataFrame(columns = columnNames, index  = indices)

        for corr in proteinCorrelations:
            prot1 = corr.prot1
            prot2 = corr.prot2
            #groundtruth
            dataTruth = self.groundtruth.copy()
            dataTruth = dataTruth[['sample_id', 'ab_id', 'ab_count_normalized']]
            spvalue = stats.spearmanr(dataTruth.loc[dataTruth.ab_id == prot1, ["ab_count_normalized"]], dataTruth.loc[dataTruth.ab_id == prot2, ["ab_count_normalized"]]).correlation
            correlations.loc[prot1 + '_' + prot2, "Groundtruth"] = float(spvalue)

            for key in self.data:
                if("Groundtruth" in key):
                    continue
                data = self.data[key].copy()
                data = data[['sample_id', 'ab_id', 'ab_count_normalized']]
                spvalue = stats.spearmanr(data.loc[data.ab_id == prot1, ["ab_count_normalized"]], data.loc[data.ab_id == prot2, ["ab_count_normalized"]]).correlation
                correlations.loc[prot1 + '_' + prot2, key] = float(spvalue)

        #correlations.astype(np.float64)
        cols = correlations.columns.values.tolist()
        correlations[cols] = correlations[cols].apply(pd.to_numeric, errors='coerce', axis=1)
        correlations.round(2)
        wantedCorrelations = correlations.copy()

        #ax= plt.subplot()
        #sns.heatmap(correlations, annot=True, annot_kws={"size": 5}, cmap="viridis", vmin = 0, vmax = 1) # font size
        #ax.set_ylabel('AB Correlations', fontsize = 10)
        #ax.set_title('AB Correlations', fontsize = 10)
        #ax.vlines([1], *ax.get_ylim(), linewidth = 2, color = "red")
        #plt.xticks(rotation = 45, ha='right', rotation_mode='anchor', fontsize = 8)
        #plt.yticks(rotation = 45, ha='right', rotation_mode='anchor', fontsize = 8)
        #plt.legend( fontsize='8', title_fontsize = '10') # for legend text
        #plt.tight_layout()
        #plt.savefig(self.folder_path + "Overview/WantedCorrelationHeatmap.png", dpi=199)
        #plt.close()

        #UNWANTED CORRELATIONS
        #for the proteinCorrelations calcualte the spearman correlation for groundtruth and norm
        unwantedProteins = self.groundtruth["ab_id"].unique()
        unwantedProteinList = [x for x in unwantedProteins if x not in wantedVarProteinList]

        negProteinCorrelations = []
        for prot1, prot2 in itertools.combinations(unwantedProteinList, 2):
            cor = (prot1, prot2)
            negProteinCorrelations.append(cor)

        sampleNumCorrelations = min(50, len(negProteinCorrelations))
        negProteinCorrelations = random.sample(negProteinCorrelations, sampleNumCorrelations)
        indices = []
        columnNames = ["Groundtruth"]
        for corr in negProteinCorrelations:
            prot1 = corr[0]
            prot2 = corr[1]
            indices.append(prot1 + '_' + prot2)
        for key in self.data:
            if("Groundtruth" in key):
                continue
            columnNames.append(key)
        correlations = pd.DataFrame(columns = columnNames, index  = indices)

        for corr in negProteinCorrelations:
            prot1 = corr[0]
            prot2 = corr[1]
            #groundtruth
            dataTruth = self.groundtruth.copy()
            dataTruth = dataTruth[['sample_id', 'ab_id', 'ab_count_normalized']]
            spvalue = stats.spearmanr(dataTruth.loc[dataTruth.ab_id == prot1, ["ab_count_normalized"]], dataTruth.loc[dataTruth.ab_id == prot2, ["ab_count_normalized"]]).correlation
            correlations.loc[prot1 + '_' + prot2, "Groundtruth"] = float(spvalue)

            for key in self.data:
                if("Groundtruth" in key):
                    continue
                data = self.data[key].copy()
                data = data[['sample_id', 'ab_id', 'ab_count_normalized']]
                spvalue = stats.spearmanr(data.loc[data.ab_id == prot1, ["ab_count_normalized"]], data.loc[data.ab_id == prot2, ["ab_count_normalized"]]).correlation
                correlations.loc[prot1 + '_' + prot2, key] = float(spvalue)

        #correlations.astype(np.float64)
        cols = correlations.columns.values.tolist()
        correlations[cols] = correlations[cols].apply(pd.to_numeric, errors='coerce', axis=1)
        correlations.round(2)

        correlations = pd.concat([wantedCorrelations, correlations])

        ax= plt.subplot()
        sns.heatmap(correlations, annot=False, annot_kws={"size": 5}, cmap="bwr", vmin = -1, vmax = 1, cbar_kws={"shrink": .70}) # font size
        # labels, title and ticks
        ax.set_ylabel('AB Correlations', fontsize = 10)
        ax.set_title('AB Correlations', fontsize = 10)
        ax.vlines([1], *ax.get_ylim(), linewidth = 2, color = "red")
        ax.hlines([len(wantedCorrelations)], *ax.get_xlim(), linewidth = 2, color = "red")
        plt.xticks(rotation = 45, ha='right', rotation_mode='anchor', fontsize = 8)
        plt.yticks(rotation = 45, ha='right', rotation_mode='anchor', fontsize = 6)
        plt.tight_layout()
        plt.savefig(self.folder_path + "Overview/UnwantedCorrelationHeatmap.png", dpi=199)
        plt.close()

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
            data = self.data[key].copy()
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

    def validate_treatment_effect(self, diffExProteins, diffExFactors):

        abIdConditionToFactorMap = self.__parseTreatmentConditionInfluence(diffExProteins, diffExFactors)
        treatmentProts = []
        #get all diff proteins
        for vec in diffExProteins:
            treatmentProts.extend(vec)
        
        fig, axs = plt.subplots(len(self.data), len(treatmentProts),figsize=(17,12))
        i = 0
        for key in self.data:
            data = self.data[key].copy()
            data = data[data.ab_id.isin(treatmentProts)]
            j = 0
            for prot in treatmentProts:
                dict = {}
                dataTmp = data[data["ab_id"] == prot]
                for condition in data.cluster_id.unique():
                    condition = str(condition)
                    dict[condition] = data.loc[data["cluster_id"] == condition, "ab_count_normalized"]
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




