from sklearn.neighbors import kneighbors_graph
import numpy as np
import pandas as pd

#build knn graph
def build_knn_graph(data, knnSize):

    #pivot data into numpy array of SAMPLE * FEATURES
    pivotData = data.pivot(index='sample_id', columns='ab_id', values='ab_count_normalized')
    # Convert the pivot DataFrame to a NumPy array
    numpyArray = pivotData.to_numpy()
    graph = kneighbors_graph(numpyArray, knnSize, mode='connectivity', include_self=True)

    return(graph.toarray())

#calcualte the KNN overlap between two graphs
#use manhattan distance
#input data has following important columns: ab_id, sample_id, ab_count_normalized (also Groundtruth)

#NUMPY ARRAY: [[a,b],[c,d]]
#sum(axis = 1) per row ([a,b] is a row)
def calculate_knn_overlap(normData, groudTruthData, knnSize=20):
    #graph is in a form: row is sample: and set values if a column is among neighbors
    #=> rows sum to knnSize, columns not necessary
    normGraph = build_knn_graph(normData, knnSize)
    trueGraph = build_knn_graph(groudTruthData, knnSize)

    #calcualte overlap
    #matrix element wise multiply (get same set ones: same enighbors) // substract row number for diagonal ones // divide by KNN number and divide by cells (rows)
    overlapMean = ( np.sum(np.multiply(normGraph, trueGraph)) - np.shape(trueGraph)[0])/ (knnSize * np.shape(trueGraph)[0])
    overlapDist = np.sum(np.multiply(normGraph, trueGraph), axis = 1) / knnSize

    pivotDataG = groudTruthData.pivot(index='sample_id', columns='ab_id', values='ab_count_normalized')
    pivotDataN = groudTruthData.pivot(index='sample_id', columns='ab_id', values='ab_count_normalized')

    assert (pivotDataG.index == pivotDataN.index).all(), "order of samples is not the same in KNN Overlap calculation in KnnSimilarity.py"
    overlapDistDf = pd.DataFrame({"sample_id":pivotDataG.index, "KnnOverlap":overlapDist})

    return(overlapDistDf)
    