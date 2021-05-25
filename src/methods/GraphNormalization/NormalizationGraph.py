import networkx as nx
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statistics
from itertools import combinations
from itertools import combinations
from numpy import mean
from numpy import var
from math import sqrt

# function to calculate Cohen's d for independent samples
def cohend(d1, d2):
	# calculate the size of samples
	n1, n2 = len(d1), len(d2)
	# calculate the variance of the samples
	s1, s2 = var(d1, ddof=1), var(d2, ddof=1)
	# calculate the pooled standard deviation
	s = sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
	# calculate the means of the samples
	u1, u2 = mean(d1), mean(d2)
	# calculate the effect size
	return (u1 - u2) / s

class NormalizationGraph:  

    #data has following structure:
    #<ab_id | ab_count | batch_id | sample_id>
    def __build_graph_from_data(self, corr_method = "spearman"):

        #generate a list of sp correlations
        data_table = self.data.loc[:,['sample_id','ab_id', 'ab_count']]

        corr_table = pd.DataFrame()
        if(corr_method == "spearman"):
            data_table=data_table.pivot(index = "ab_id", columns='sample_id', values='ab_count')
            sp = stats.spearmanr(data_table, axis = 1)
            corr_table = pd.DataFrame(sp.correlation, columns = data_table.index, index = data_table.index)
            corr_table.reset_index(inplace=True)
            corr_table = (corr_table.melt(id_vars=['ab_id'],  var_name='ab_id_2', value_name='correlation'))
        elif(corr_method == "pearson"):
            data_table = data_table.pivot(index = "sample_id", columns='ab_id', values='ab_count')
            corr_table = pd.DataFrame(columns = ["ab_id", "ab_id_2", "correlation"])
            for x,y in (combinations(data_table.columns,2)):
                corr_table = corr_table.append(pd.DataFrame({"ab_id":[x], "ab_id_2":[y], "correlation":[stats.pearsonr(data_table[x], data_table[y])[0]]}))

        #filter duplicates from correlation matrix
        filter_cols = corr_table.filter(like='ab_id').values
        filtered_indices = pd.DataFrame(np.sort(filter_cols, axis=1)).duplicated()
        corr_table = corr_table[~filtered_indices]
        corr_table=corr_table[corr_table['ab_id'] != corr_table['ab_id_2']]

        #make a list of all edges
        edge_list = list()
        for index, row in corr_table.iterrows():
            if(row['correlation'] > 0.7):
                edge_list.append((row['ab_id'], row['ab_id_2'], {'weight': row['correlation']}))

        #conatruct graph from edge list
        G = nx.Graph()
        G.add_edges_from(edge_list)

        #nx.draw(G, with_labels=True, font_weight='bold')
        #plt.show()

        return(G)

    def __normalize_by_library_size(self, data):
        data["ab_count_compositional"] = data.apply(lambda row : (row["ab_count"]/sum(data.loc[data["sample_id"]==row["sample_id"],"ab_count"])), axis = 1)
        return(data)
    def __calculate_library_size(self, data):
        data["lib_size"] = data.apply(lambda row : (sum(data.loc[data["sample_id"]==row["sample_id"],"ab_count"])), axis = 1)
        return(data)
    #initialize graph for normalization
    #vertexes are protein abundancies, edges their correlations
    def __init__(self, data):

        self.data = self.__normalize_by_library_size(data)
        self.data=data
        #Graph
        self.G = self.__build_graph_from_data()

    def list_of_no_correlated_samples(self, clique_list):
        new_clique_list = list()
        for clique in clique_list:
            cluster_data_table = self.data.loc[:,['sample_id','ab_id', 'ab_count_compositional', 'cluster_id']]
            treatments = cluster_data_table["cluster_id"].unique()
            treatment_dict = dict()
            abs = clique.copy()
            test_dict = dict()
            for t in treatments:
                treatment_table = cluster_data_table[cluster_data_table["cluster_id"] == t]
                treatment_table = treatment_table.pivot(index = "sample_id", columns='ab_id', values='ab_count_compositional')
                test_dict[t] = treatment_table
            for ab in abs:
                for x,y in (combinations(treatments,2)):
                    ttest, pval = stats.ttest_ind(test_dict[x][ab], test_dict[y][ab])
                    cohen_d = cohend(test_dict[x][ab], test_dict[y][ab])
                    if(pval<0.05 and cohen_d >= 0.5):
                        clique.remove(ab)
                        break
            new_clique_list.append(clique)
        return(new_clique_list)  

    #return a table with ID and norm_score column
    def get_normalized_score(self):
        #find max_clique
        clique_list = list(nx.find_cliques(self.G))
        max_clique = list()
        max = 0

        clique_list = self.list_of_no_correlated_samples(clique_list)

        for clique in clique_list:
            if(len(clique) > max):
                max = len(clique)
                max_clique = clique
        print(max)
        print("Found Clique: ")
        print(max_clique)

        #normalize by using these feaures
        feature_mean = dict()
        for feature in max_clique:
            ab_values = self.data[self.data["ab_id"] == feature]
            mean = statistics.mean(ab_values["ab_count"])
            feature_mean[feature] = mean

        normalized_data_frame = pd.DataFrame()
        for sample in np.unique(self.data["sample_id"]):
            sample_data = self.data[self.data["sample_id"] == sample].copy()
            avg_scaling_factor = 0
            for feature in feature_mean:
                mean_value = feature_mean[feature]
                sample_value = sample_data[sample_data["ab_id"] == feature]
                avg_scaling_factor += (mean_value/ sample_value.iloc[0]["ab_count"])

            avg_scaling_factor = avg_scaling_factor/len(feature_mean)
            sample_data["ab_count_normalized"] = np.log((sample_data["ab_count"]*avg_scaling_factor))
            normalized_data_frame = normalized_data_frame.append(sample_data)
        return(normalized_data_frame)

    #return a table with ID and norm_score column
    #normalize over all clusters ith more than 20 nodes and average over those
    def get_normalized_score_multiClique(self):
        #find max_clique
        clique_list = list(nx.find_cliques(self.G))
        clique_list = self.list_of_no_correlated_samples(clique_list)

        count = 0
        sorted_clique = sorted(clique_list, key=len, reverse=True)
        new = list()
        for el in sorted_clique:
            if(len(el) > 7):
                new.append(el)
        clique_list = new
        print(len(clique_list))

        weight_list = list()
        normVector_list = list()
        for clique in clique_list:
            normalized_data_frame = pd.DataFrame()
            print(str(count) + " " + str(len(clique)))
            #normalize by using these feaures
            feature_mean = dict()
            for feature in clique:
                ab_values = self.data[self.data["ab_id"] == feature]
                mean = statistics.mean(ab_values["ab_count"])
                feature_mean[feature] = mean

            sample_vector = np.array([])
            for sample in np.unique(self.data["sample_id"]):
                sample_data = self.data[self.data["sample_id"] == sample].copy()
                avg_scaling_factor = 0
                for feature in feature_mean:
                    mean_value = feature_mean[feature]
                    sample_value = sample_data[sample_data["ab_id"] == feature]
                    avg_scaling_factor += (mean_value/ sample_value.iloc[0]["ab_count"])

                avg_scaling_factor = avg_scaling_factor/len(feature_mean)
                sample_vector = np.concatenate((sample_vector, np.log((sample_data["ab_count"]*avg_scaling_factor))))
                normalized_data_frame = normalized_data_frame.append(sample_data)
            
            normVector_list.append(sample_vector)
            weight_list.append(len(clique))      
            count += 1
        weight_list = np.array(weight_list)
        normVector_list = np.stack(normVector_list)

        norm_matrix = normVector_list.T*weight_list
        norm_matrix_t = sum(norm_matrix.T)

        normalized_data_frame["ab_count_normalized"] = norm_matrix_t/sum(weight_list)
        return(normalized_data_frame)    

