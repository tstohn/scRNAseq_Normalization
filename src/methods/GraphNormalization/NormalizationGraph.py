import networkx as nx
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class NormalizationGraph:

    #data has following structure:
    #<ab_id | ab_count | batch_id | sample_id>
    def build_graph_from_data(self):

        #generate a list of sp correlations
        data_table = self.data.loc[:,['sample_id','ab_id', 'ab_count']]

        data_table=data_table.pivot(index = "ab_id", columns='sample_id', values='ab_count')

        sp = stats.spearmanr(data_table, axis = 1)
        corr_table = pd.DataFrame(sp.correlation, columns = data_table.index, index = data_table.index)
        corr_table.reset_index(inplace=True)
        corr_table = (corr_table.melt(id_vars=['ab_id'],  var_name='ab_id_2', value_name='pearson'))

        #filter duplicates from correlation matrix
        filter_cols = corr_table.filter(like='ab_id').values
        filtered_indices = pd.DataFrame(np.sort(filter_cols, axis=1)).duplicated()
        corr_table = corr_table[~filtered_indices]
        corr_table=corr_table[corr_table['ab_id'] != corr_table['ab_id_2']]

        edge_list = list()
        for index, row in corr_table.iterrows():
            if(row['pearson'] > 0.9):
                edge_list.append((row['ab_id'], row['ab_id_2'], {'weight': row['pearson']}))

        print(edge_list[0])
        G = nx.Graph()
        G.add_edges_from(edge_list)

        nx.draw(G, with_labels=False, font_weight='bold')
        plt.savefig("graph.png")

        clique_list = list(nx.find_cliques(G))
        max_clique = list()
        max = 0
        for clique in clique_list:
            if(len(clique) > max):
                max = len(clique)
                max_clique = clique

        print(max)
        print(max_clique)        

        return(G)

    #initialize graph for normalization
    #vertexes are protein abundancies, edges their correlations
    def __init__(self, data):

        self.data = data
        #Graph
        self.G = self.build_graph_from_data()
 
    #return a table with ID and norm_score column
    def get_normalized_score(self):
        #find max_clique

        #check clique for uniform distribution

        #normalize by using these feaures
        return(1)

