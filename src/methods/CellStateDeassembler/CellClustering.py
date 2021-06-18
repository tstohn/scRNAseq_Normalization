import numpy as np
import meld

class Cluster:
    def __init__(self, data):
        self.data = data

#represent cells by a function of cell state (cell stae is a continuous space)
#model this space as a manifold, create a manufold for each treatment and make a mapping of 
#a manifold distriution function describing all possible cell states to this function for a different treatment

#put those mappings into a auto-encoder to elarn the hidden continuous space of cell state points for different treatments
class CellClustering:
    def __init__(self, data):
        self.data = data

    def meld_algorithm(self):
        pre_data = self.data.filter(items=['sample_id', 'ab_count_normalized', 'cluster_id', 'ab_id'])
        pre_data = pre_data.filter(items=['sample_id', 'ab_count_normalized', 'ab_id', 'cluster_id'])
        pre_data = pre_data.pivot(index = 'sample_id', columns = 'ab_id', values = 'ab_count_normalized')
        data = pre_data.sort_values(by = ['sample_id']).values

        pre_labels = self.data.filter(items = ['cluster_id', 'sample_id'])
        pre_labels = pre_labels.set_index('sample_id')
        sample_labels = pre_labels.sort_values(by = ['sample_id']).values

        sample_densities = meld.MELD().fit_transform(data, sample_labels)
        #do for all cluster_ids
        egf_likelihoods = meld.utils.normalize_densities(sample_densities)["EGF"]
        iRSK_EGF_likelihoods = meld.utils.normalize_densities(sample_densities)["iRSK_EGF"]
        ip70S6K_EGF_likelihoods = meld.utils.normalize_densities(sample_densities)["ip70S6K_EGF"]
