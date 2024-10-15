from sklearn.neighbors import kneighbors_graph
import numpy as np
import pandas as pd
from scipy.stats import zscore

#build knn graph
def build_knn_graph(data, knnOverlap, knnMetric):

    #pivot data into numpy array of SAMPLE * FEATURES
    pivotData = data.pivot(index='sample_id', columns='ab_id', values='ab_count_normalized')
    pivotData = pivotData.apply(zscore)

    # Convert the pivot DataFrame to a NumPy array
    numpyArray = pivotData.to_numpy()

    #distances: cityblock, cosine, euclidean (re-write metric for 'strange' sklearn name)
    if(knnMetric == "manhattan"):
        knnMetricTmp = "cityblock"
    else:
        knnMetricTmp = knnMetric
    graph = kneighbors_graph(numpyArray, knnOverlap, mode='connectivity', include_self=True, metric=knnMetricTmp)

    return(graph.toarray())

#calcualte the KNN overlap between two graphs
#use manhattan distance
#input data has following important columns: ab_id, sample_id, ab_count_normalized (also Groundtruth)

#NUMPY ARRAY: [[a,b],[c,d]]
#sum(axis = 1) per row ([a,b] is a row)
def calculate_knn_overlap(normData, groudTruthData, knnOverlap, knnMetric):
    #graph is in a form: row is sample: and set values if a column is among neighbors
    #=> rows sum to knnSize, columns not necessary
    
    normGraph = build_knn_graph(normData, knnOverlap, knnMetric)
    trueGraph = build_knn_graph(groudTruthData, knnOverlap, knnMetric)

    #calcualte overlap
    #matrix element wise multiply (get same set ones: same enighbors) // substract row number for diagonal ones // divide by KNN number and divide by cells (rows)
    overlapMean = ( np.sum(np.multiply(normGraph, trueGraph)) - np.shape(trueGraph)[0])/ (knnOverlap * np.shape(trueGraph)[0])
    overlapDist = np.sum(np.multiply(normGraph, trueGraph), axis = 1) / knnOverlap

    pivotDataG = groudTruthData.pivot(index='sample_id', columns='ab_id', values='ab_count_normalized')
    pivotDataN = groudTruthData.pivot(index='sample_id', columns='ab_id', values='ab_count_normalized')

    assert (pivotDataG.index == pivotDataN.index).all(), "order of samples is not the same in KNN Overlap calculation in KnnSimilarity.py"
    overlapDistDf = pd.DataFrame({"sample_id":pivotDataG.index, "KnnOverlap":overlapDist, "KnnThreshold":knnOverlap})

    return(overlapDistDf)
    