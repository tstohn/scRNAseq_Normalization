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
import scanpy as sc

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
    def __build_graph_from_data(self, corr_method, corr_threshold):

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
            if(abs(row['correlation']) > corr_threshold):
                edge_list.append((row['ab_id'], row['ab_id_2'], {'weight': row['correlation']}))

        #conatruct graph from edge list
        G = nx.Graph()
        G.add_edges_from(edge_list)

        #nx.draw(G, with_labels=True, font_weight='bold')
        #plt.show()

        return(G)

    def normalize_by_library_size(self):
        summed_data=self.data.groupby(['sample_id'])['ab_count_normalized'].sum().reset_index()
        summed_data.rename(columns = {'ab_count_normalized': 'ab_count_compositional'}, inplace = True)
        summed_data['sample_id'].astype(str)
        self.data['sample_id'].astype(str)
        self.data = self.data.merge(summed_data)
        self.data["ab_count_compositional"] = self.data.ab_count_normalized / self.data.ab_count_compositional
        return(self.data)

    def __calculate_library_size(self, data):
        data["lib_size"] = data.apply(lambda row : (sum(data.loc[data["sample_id"]==row["sample_id"],"ab_count"])), axis = 1)
        return(data)
    #initialize graph for normalization
    #vertexes are protein abundancies, edges their correlations
    def __init__(self, data, corr_method = "spearman", corr_threshold=0.7):

        self.data=data
        self.data["ab_count_normalized"] = self.data["ab_count"]
        #Graph

        self.G = self.__build_graph_from_data(corr_method, corr_threshold)

    def list_of_no_correlated_samples(self, clique_list, p_val, cohend_val):
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
                    print(treatments)
                    print(str(x) + "_" + str(y) + ":" + str(abs))
                    ttest, pval = stats.ttest_ind(test_dict[x][ab], test_dict[y][ab])
                    cohen_d = cohend(test_dict[x][ab], test_dict[y][ab])
                    if(pval<p_val and cohen_d >= cohend_val):
                        clique.remove(ab)
                        break
            new_clique_list.append(clique)
        return(new_clique_list)  

    #return a table with ID and norm_score column
    def get_normalized_score(self, p_val=0.05, cohend_val=0.5, take_log=False):
        #find max_clique
        clique_list = list(nx.find_cliques(self.G))
        max_clique = list()
        max = 0
        clique_list = self.list_of_no_correlated_samples(clique_list, p_val, cohend_val)

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
            mean = statistics.mean(ab_values["ab_count_normalized"])
            feature_mean[feature] = mean

        normalized_data_frame = pd.DataFrame()
        for sample in np.unique(self.data["sample_id"]):
            sample_data = self.data[self.data["sample_id"] == sample].copy()
            avg_scaling_factor = 0
            zero_scalings = 0
            for feature in feature_mean:
                mean_value = feature_mean[feature]
                sample_value = sample_data[sample_data["ab_id"] == feature]
                if(sample_value.iloc[0]["ab_count_normalized"] == 0 or mean_value == 0):
                    zero_scalings+=1
                else:
                    avg_scaling_factor += (mean_value/ sample_value.iloc[0]["ab_count_normalized"])

            assert(len(feature_mean) != zero_scalings)
            avg_scaling_factor = avg_scaling_factor/(len(feature_mean)-zero_scalings)
            if(take_log):
                sample_data["ab_count_normalized"] = np.log(sample_data["ab_count_normalized"]*avg_scaling_factor)
            else:
                sample_data["ab_count_normalized"] = (sample_data["ab_count_normalized"]*avg_scaling_factor)
            normalized_data_frame = normalized_data_frame.append(sample_data)

        self.data = normalized_data_frame
        return(normalized_data_frame)

    #return a table with ID and norm_score column
    #normalize over all clusters ith more than 20 nodes and average over those
    #NOT USED at the moment
    def get_normalized_score_multiClique(self, take_log=False):
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
                zero_scalings = 0
                for feature in feature_mean:
                    mean_value = feature_mean[feature]
                    sample_value = sample_data[sample_data["ab_id"] == feature]
                    if(sample_value.iloc[0]["ab_count"] == 0 or mean_value == 0):
                        zero_scalings+=1
                    else:
                        avg_scaling_factor += (mean_value/ sample_value.iloc[0]["ab_count"])

                assert(len(feature_mean) != zero_scalings)
                avg_scaling_factor = avg_scaling_factor/(len(feature_mean)-zero_scalings)
                if(take_log):
                    sample_vector = np.concatenate(sample_vector, np.log(sample_data["ab_count"]*avg_scaling_factor))
                else:
                    sample_vector = np.concatenate((sample_vector, (sample_data["ab_count"]*avg_scaling_factor)))
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

    def remove_batch_effect(self):

        #pivot data (samples*ABs)
        combat_data = self.data[["sample_id", "ab_id", "batch_id", "ab_count_normalized"]]
        combat_data = combat_data.pivot(index=["sample_id","batch_id"], columns="ab_id", values="ab_count_normalized")

        #generate an AnnData matrix (index is batch and sample id, 
        #batch is then taken as column for observation matrix as well)
        index_array = [np.array(combat_data.index.get_level_values(0)),
                       np.array(combat_data.index.get_level_values(1))]
        var_df = pd.DataFrame(index=combat_data.columns.values)
        obs_df = pd.DataFrame(index=pd.MultiIndex.from_arrays(index_array))
        obs_df["batch_id"] = combat_data.index.get_level_values(1)
        obs_df["batch_id"] = obs_df["batch_id"].astype("category")
        anndata_matrix = sc.AnnData(X=combat_data, var = var_df, obs = obs_df)

        #call scanpy
        sc.pp.combat(anndata_matrix, key="batch_id", inplace = True)

        #revert AnnData matrix back to a dataframe and merge with origional one
        reverted_ann_matrix = pd.DataFrame(anndata_matrix.X, index = anndata_matrix.obs.index, columns = anndata_matrix.var.index)
        reindexed_matrix = reverted_ann_matrix.reset_index()
        unpivoted_matrix = reindexed_matrix.melt(id_vars=["level_0", "level_1"], var_name="ab_id")
        normalized_matrix = unpivoted_matrix.rename(columns = {'level_0':'sample_id', 'level_1':'batch_id', 'value':'ab_count_normalized'})
        self.data.drop('ab_count_normalized', inplace = True, axis=1)
        result = self.data.merge(normalized_matrix)

        self.data = result
        return(result)