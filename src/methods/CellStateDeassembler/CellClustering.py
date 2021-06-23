from numba.cuda.api import per_thread_default_stream
import numpy as np
import pandas as pd
import meld
import phate
import scprep
import matplotlib.pyplot as plt
from sklearn import cluster
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
import umap
from sklearn.decomposition import PCA

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

    def z_score(self, df):
        # copy the dataframe
        df_std = df.copy()
        # apply the z-score method
        for column in df_std.columns:
            df_std[column] = (df_std[column] - df_std[column].mean()) / df_std[column].std()
            
        return df_std

    def umap_function(self, data):
        standard_embedding = umap.UMAP(random_state=42).fit_transform(data)
        kmeans_labels = cluster.KMeans(n_clusters=10).fit_predict(standard_embedding)
        return kmeans_labels

    def mean_shift_function(self, data):
        print(data)
        bandwidth = cluster.estimate_bandwidth(data, quantile=0.05)
        #bandwidth = 6.5
        clustering = cluster.MeanShift(bandwidth=bandwidth).fit(data)
        all_labels = clustering.labels_
        return all_labels

    def kde_loglikely_function(self, data):
        kde = KernelDensity(bandwidth=1, kernel='epanechnikov')
        kde.fit(data)
        # score_samples returns the log of the probability density
        logprob = kde.score_samples(data)
        return logprob

    def melt_egf_function(self, data):

        pre_labels = self.data.filter(items = ['cluster_id', 'sample_id'])
        pre_labels = pre_labels.drop_duplicates()
        pre_labels = pre_labels.set_index('sample_id')
        sample_labels = pre_labels.sort_values(by = ['sample_id'])
        sample_labels = sample_labels.values

        sample_densities = meld.MELD().fit_transform(data, sample_labels)
        #do for all cluster_ids
        egf_likelihoods = meld.utils.normalize_densities(sample_densities)["EGF"]
        
        return egf_likelihoods.values

    def get_kmean_labels_for_all(self, data, label_func, zscore = True):
        all_labels = None

        data = data.filter(items=['sample_id', 'ab_count_normalized', 'cluster_id', 'ab_id'])
        data = data.pivot(index = 'sample_id', columns = 'ab_id', values = 'ab_count_normalized')
        data = data.sort_values(by = ['sample_id'])
        if(zscore):
            data = self.z_score(data)

        all_labels = label_func(data.values)
        data["all_labels"] = all_labels
        sample_to_label_mapping = data.filter(items=['sample_id', 'all_labels'])
        
        return sample_to_label_mapping

    #writ epivotted data with sample-label mapping to tsv file
    def write_data_frame(self, data, labels):
        data = data.reset_index()
        data = pd.melt(data, id_vars=['sample_id'])
        out_data = pd.merge(data, labels, on=["sample_id"])
        out_data.to_csv("/Users/t.stohn/Desktop/TMP/OUTDATA.tsv")

    #for now: calcualte mean shift clusters of different treatments
    #OUTLOOK: map clusters of different treatments
    def cluster(self):
        cluster_ids = self.data.cluster_id.unique()
        #cluster_ids = np.append(cluster_ids, "ALL")

        cluster_ids = np.insert(cluster_ids, 0, "ALL")
        mean_shift_data = self.data.filter(items=['sample_id', 'ab_count_normalized', 'cluster_id', 'ab_id'])

        #call function to generte labels for kmean over all samples
        labels = self.get_kmean_labels_for_all(mean_shift_data, self.mean_shift_function)
        mean_shift_data = pd.merge(mean_shift_data, labels, on=["sample_id"])
        print(len(np.unique(labels)))
        print(labels.value_counts()
)

        for single_cluster in cluster_ids:

            if(not single_cluster=="ALL"):
                pre_data = mean_shift_data.loc[mean_shift_data["cluster_id"] == single_cluster]
            else:
                pre_data = mean_shift_data

            pre_data = pre_data.pivot(index = ['sample_id', 'all_labels', 'cluster_id'] , columns = 'ab_id', values = 'ab_count_normalized')
            pre_data = pre_data.reset_index("all_labels")
            pre_data = pre_data.reset_index("cluster_id")
            pre_data = pre_data.sort_index(ascending = True)

            all_labels = pre_data["all_labels"]
            pre_data = pre_data.drop(columns = 'all_labels')

            cluster_labels = pre_data["cluster_id"]
            pre_data = pre_data.drop(columns = 'cluster_id')
            #pre_data = self.z_score(pre_data)

            if(single_cluster=="ALL"):
                self.write_data_frame(pre_data, labels)

            data = pre_data.sort_values(by = ['sample_id']).values
            all_labels = all_labels.values
            cluster_labels = cluster_labels.values
            cluster_labels = pd.factorize(cluster_labels)[0].tolist()

            #cluster subset on MEAN SHIFT
            bandwidth = cluster.estimate_bandwidth(data, quantile=0.01, n_samples=544)
            #bandwidth = 0.68
            clustering = cluster.MeanShift(bandwidth=bandwidth).fit(data)

            #cluster subset on KMEANS
            standard_embedding = umap.UMAP(random_state=42).fit_transform(data)
            kmeans_labels = cluster.KMeans(n_clusters=10).fit_predict(standard_embedding)

            plt.figure(figsize=(16,10))
            plt.scatter(standard_embedding[:, 0], standard_embedding[:, 1], c=all_labels, s=10, cmap='Spectral');
            plt.savefig("/Users/t.stohn/Desktop/TMP/UMAP" + single_cluster + ".png")
            plt.close()

        return None

    def cluster_pca_embedding(self):
        cluster_ids = self.data.cluster_id.unique()
        #cluster_ids = np.append(cluster_ids, "ALL")

        cluster_ids = np.insert(cluster_ids, 0, "ALL")
        mean_shift_data = self.data.filter(items=['sample_id', 'ab_count_normalized', 'cluster_id', 'ab_id'])
          
        pre_data = mean_shift_data.pivot(index = ['sample_id', 'cluster_id'] , columns = 'ab_id', values = 'ab_count_normalized')
        pre_data = pre_data.reset_index("cluster_id")
        pre_data = pre_data.sort_index(ascending = True)
        cluster_labels = pre_data["cluster_id"]
        pre_data = pre_data.drop(columns = 'cluster_id')
        pre_data = self.z_score(pre_data)

        #getting the embedding of data in umap
        pca = PCA(n_components=3)
        pc_output = pca.fit_transform(pre_data)
        col_vec = []
        for i in range(0,3):
            col_vec = np.append(col_vec,("PC_"+str(i)))

        pc_data = pd.DataFrame(data = pc_output, columns = col_vec)
        pc_data['sample_id'] = pre_data.index

        data = pd.melt(pc_data, id_vars=['sample_id'])
        data = data.rename(columns={"variable": "ab_id", "value": "ab_count_normalized"})

        #data = pd.merge(data, self.data.filter(items=['sample_id', 'cluster_id']), on=["sample_id"])

        #call function to generte labels for kmean over all samples
        print(data)
        labels = self.get_kmean_labels_for_all(data, self.mean_shift_function, zscore=False)

        #write data to tsv
        out_data = pd.merge(mean_shift_data, labels, on=["sample_id"])
        out_data.to_csv("/Users/t.stohn/Desktop/TMP/PCA_LABELLED_OUTPUT.tsv")

        return None

    #Calculate KDE for the untreated(baseline) population
    #calculate difference for all treated cells compared to baseline density
    #check which cells have diverged
    def cluster_KDE_diff(self):
        #generate KDE of baseline data (in this example case EGF)
        cluster = "ALL"
        if(cluster == "EGF"):
            data = self.data.loc[self.data["cluster_id"] == cluster]
        else:
            data = self.data
        
        data = data.filter(items=['sample_id', 'ab_count_normalized', 'cluster_id', 'ab_id'])
        data = data.pivot(index = 'sample_id', columns = 'ab_id', values = 'ab_count_normalized')
        data = data.sort_values(by = ['sample_id'])
        X = self.z_score(data)
        X = X.values

        #estimate a good bandwidth
        params = {'bandwidth': np.linspace(-1, 1, 20)}
        grid = GridSearchCV(KernelDensity(), params, cv=20)
        grid.fit(X)

        print("best bandwidth: {0}".format(grid.best_estimator_.bandwidth))

        #calculate KDE
        kde = grid.best_estimator_
        estimate = kde.score_samples(X)

        #instantiate and fit the KDE model
        kde = KernelDensity(bandwidth=6.5, kernel='gaussian')
        kde.fit(X)
        # score_samples returns the log of the probability density
        logprob = kde.score_samples(X)
        print(logprob)

        return logprob

        #evaluate accuracy of estimator
