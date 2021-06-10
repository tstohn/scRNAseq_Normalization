import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from random import randint
from scipy import stats

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
        if not os.path.exists("bin/BENCHMARKED_DATASETS/"+folder_name + "/ConfusionMatrices"):
            os.mkdir("bin/BENCHMARKED_DATASETS/"+folder_name + "/ConfusionMatrices")
        if not os.path.exists("bin/BENCHMARKED_DATASETS/"+folder_name + "/MMDMatrices"):
            os.mkdir("bin/BENCHMARKED_DATASETS/"+folder_name + "/MMDMatrices")
        self.results = open("bin/BENCHMARKED_DATASETS/"+folder_name+"/results.tsv", "w+")
        self.results.write("NORMALIZATION_METHOD" + "\t" + "CLASSIFICATION_METHOD" + "\t" + "ACCURACY_MEAN" + "\t" + "ACCURACY_SD" + "\n")

        self.sp_results = open("bin/BENCHMARKED_DATASETS/"+folder_name+"/results_spearman.tsv", "w+")
        self.sp_results.write("NORMALIZATION_METHOD" + "\t" + "SPEARMAN_CORRELATION_MEAN" + "\t" + "SPEARMAN_PVALUE_MEAN" + "\n")

        self.dataset_name = folder_name
        self.folder_path = ("bin/BENCHMARKED_DATASETS/"+folder_name+"/")

    #def __del__(self):
        #self.results.close()

    def __calculate_scores(self, model, params, X, y, global_scores, lock, random_factor=None):
        cv_inner = KFold(n_splits=10, shuffle=True, random_state=random_factor)
        search = RandomizedSearchCV(model, params, n_iter = 30, scoring='accuracy', n_jobs=1, cv=cv_inner, refit=True, random_state=random_factor)
        cv_outer = KFold(n_splits=5, shuffle=True, random_state=random_factor)
        scores = cross_val_score(search, X, y, scoring='accuracy', cv=cv_outer, n_jobs=-1)
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
        scores = cross_val_score(search, X, y, scoring='accuracy', cv=cv_outer, n_jobs=-1)
        #plot a confusion matrix
        y_pred = cross_val_predict(search, X, y, cv=cv_outer, n_jobs=-1)
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
        plt.savefig(self.folder_path + "ConfusionMatrices/" + data_name + method_string + "_CFMatrix.png")
        plt.close()

        print('%s => %s Accuracy[%s] : %.3f (%.3f)' % (self.dataset_name, method_string, data_name, np.mean(global_scores), np.std(global_scores)))
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

    def calculate_wanted_and_unwanted_variance(self):
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
            plt.savefig(self.folder_path + "MMDMatrices/" + key + ".png")
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
        plt.savefig(self.folder_path + data_name + "_output.png")
        plt.close()
        
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

    def ab_spearman_correlation(self, filter = ""):
        plt.figure(figsize=(16,10))
        lengend_labels = []
        for key in self.data:
            data = self.data[key]

            #speacial line for a data set where phospho proteins were excluded due to expected correlations
            if(filter != ""):
                data = data[data['ab_type'] == filter]

            data_subset = data.loc[:,['sample_id','ab_id', 'ab_count_normalized']]
            d_pivot=data_subset.pivot(index = "ab_id", columns='sample_id', values='ab_count_normalized')
            
            try:
                sp, p = stats.spearmanr(d_pivot, axis = 1)
                dim=len(sp[0])
            except:
                drop_list = []
                for i in range(len(d_pivot.index)):
                    row = d_pivot.iloc[i].values.tolist()
                    if(all(elem == row[0] for elem in row)):
                        drop_list.append(i)
                d_pivot = d_pivot.drop(d_pivot.index[drop_list], axis = 0)
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
        plt.savefig(self.folder_path +  "spearman_correlations_" + filter + ".png")
        plt.close()






