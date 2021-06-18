import numpy as np
import pandas as pd
import meld
import phate
import scprep
import matplotlib.pyplot as plt
from sklearn import cluster
import umap

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
        pre_data = pre_data.pivot(index = 'sample_id', columns = 'ab_id', values = 'ab_count_normalized')
        data = pre_data.sort_values(by = ['sample_id']).values

        pre_labels = self.data.filter(items = ['cluster_id', 'sample_id'])
        pre_labels = pre_labels.drop_duplicates()
        pre_labels = pre_labels.set_index('sample_id')
        sample_labels = pre_labels.sort_values(by = ['sample_id']).values

        sample_densities = meld.MELD().fit_transform(data, sample_labels)
        #do for all cluster_ids
        egf_likelihoods = meld.utils.normalize_densities(sample_densities)["EGF"]
        iRSK_EGF_likelihoods = meld.utils.normalize_densities(sample_densities)["iRSK_EGF"]
        ip70S6K_EGF_likelihoods = meld.utils.normalize_densities(sample_densities)["ip70S6K_EGF"]

        result_data = {'EGF': egf_likelihoods, 'iRSK_EGF': iRSK_EGF_likelihoods, 'ip70S6K_EGF': ip70S6K_EGF_likelihoods}
        result = pd.DataFrame(data=result_data)


        phateop = phate.PHATE()
        data_phate = phateop.fit_transform(data)

        #for every cluster do
        plt.figure(figsize=(16,10))
        scprep.plot.scatter2d(data_phate, c=egf_likelihoods,
                      ticks=False, figsize=(6,5),
                     vmin=0, vmax=1, cmap=meld.utils.get_meld_cmap())
        plt.savefig("/Users/t.stohn/Desktop/TMP/phate_output_EGF.png")
        plt.close()

        plt.figure(figsize=(16,10))
        scprep.plot.scatter2d(data_phate, c=iRSK_EGF_likelihoods,
                      ticks=False, figsize=(6,5),
                     vmin=0, vmax=1, cmap=meld.utils.get_meld_cmap())
        plt.savefig("/Users/t.stohn/Desktop/TMP/phate_output_iRSK.png")
        plt.close()

        plt.figure(figsize=(16,10))
        scprep.plot.scatter2d(data_phate, c=ip70S6K_EGF_likelihoods,
                      ticks=False, figsize=(6,5),
                     vmin=0, vmax=1, cmap=meld.utils.get_meld_cmap())
        plt.savefig("/Users/t.stohn/Desktop/TMP/phate_output_ip70S6K.png")
        plt.close()

        return result

    #for now: calcualte mean shift clusters of different treatments
    #OUTLOOK: map clusters of different treatments
    def mean_shift(self):
        cluster_ids = self.data.cluster_id.unique()
        cluster_ids = np.append(cluster_ids, "ALL")

        for single_cluster in cluster_ids:
            pre_data = self.data.filter(items=['sample_id', 'ab_count_normalized', 'cluster_id', 'ab_id'])

            if(not single_cluster=="ALL"):
                pre_data = pre_data.loc[pre_data["cluster_id"] == single_cluster]

            pre_data = pre_data.pivot(index = 'sample_id', columns = 'ab_id', values = 'ab_count_normalized')
            data = pre_data.sort_values(by = ['sample_id']).values

            bandwidth = cluster.estimate_bandwidth(data, quantile=0.5, n_samples=544)
            bandwidth = 0.03
            clustering = cluster.MeanShift(bandwidth=bandwidth).fit(data)
            print((clustering.labels_))

            standard_embedding = umap.UMAP(random_state=42).fit_transform(data)

            kmeans_labels = cluster.KMeans(n_clusters=10).fit_predict(standard_embedding)

            plt.figure(figsize=(16,10))
            plt.scatter(standard_embedding[:, 0], standard_embedding[:, 1], c=kmeans_labels, s=10, cmap='Spectral');
            plt.savefig("/Users/t.stohn/Desktop/TMP/UMAP" + single_cluster + ".png")
            plt.close()

        return None