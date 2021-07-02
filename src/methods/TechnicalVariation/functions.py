import sys
import pandas as pd
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy import stats

# get difference of distributions for all ABs (possible metrices: ab count and relative count)
""" return a dictionary with the average distance between all distributions of all possible tags for an ab_id:
    { <ab_id>: <average distance>} 
"""
def calculate_distribution_difference():


# rate how much a sample is an outlier to correlation
def calc_distance(p1, p2, p3):
    dist = norm(np.cross(p2-p1, p1-p3))/norm(p2-p1)
    return(dist)

def correlate_two_abs(tag1, tag2, data):
    #filter tag and select only sample and count column
    tag1_dataFrame = data[data["Ab_barcode_nr"] == tag1]
    tag1_dataFrame = tag1_dataFrame[["sample_id", "ab_count_relative"]]
    tag1_dataFrame.rename(columns = {'ab_count_relative': 'tag1'}, inplace = True)
    tag1_dataFrame["tag1"] = tag1_dataFrame["tag1"].astype(float)
    tag1_dataFrame["sample_id"] = tag1_dataFrame["sample_id"].astype(str)

    tag2_dataFrame = data[data["Ab_barcode_nr"] == tag2]
    tag2_dataFrame = tag2_dataFrame[["sample_id", "ab_count_relative"]]
    tag2_dataFrame.rename(columns = {'ab_count_relative': 'tag2'}, inplace = True)
    tag2_dataFrame["tag2"] = tag2_dataFrame["tag2"].astype(float)
    tag2_dataFrame["sample_id"] = tag2_dataFrame["sample_id"].astype(str)

    commonData = tag1_dataFrame.merge(tag2_dataFrame, on = "sample_id", how = "inner")
    tag1_data = commonData.tag1
    tag2_data = commonData.tag2

    #correlate data for tag1 and tag2
    m, b = np.polyfit(tag1_data, tag2_data, 1)
    p1 = np.array([0, b])
    p2 = np.array([1, m*1 + b])

    samples = commonData.sample_id.unique()
    columnName = str(tag1) + "_" + str(tag2)
    resultDict = {columnName : np.nan}
    resultDataFrame = pd.DataFrame(resultDict, index = samples)  

    for sample in samples:
        #the row in the dataFrame with our sample
        sample_data = commonData[commonData["sample_id"] == sample]
        #get tag1 and tag2, dataFrame has only on row (index by 0)
        p3 = np.array([sample_data.iloc[0]["tag1"], sample_data.iloc[0]["tag2"]])
        dist = calc_distance(p1, p2, p3)
        
        resultDataFrame.loc[sample, columnName] = dist

    #plt.plot(tag1_data, tag2_data, 'o')
    #plt.plot(tag1_data, m*tag1_data + b)
    #plt.show()
    #calculate distane from this correlation for every sample
    return(resultDataFrame)

""" return a matrix [sample * abType] to show the distance for every sample for all ab types to linear model"""
def distance_from_correlation(data):
    ab_types = data.ab_id.unique()
    samples = data.sample_id.unique()
    resultDataFrame = pd.DataFrame(None, index = samples)    
    for ab in ab_types:
        ab_data = data[data["ab_id"] == ab]
        tags = ab_data.Ab_barcode_nr.unique()
        #correlate all abs with each other
        result = [correlate_two_abs(tags[tag1_idx], tags[tag2_idx], data) for tag1_idx in range(len(tags)) for tag2_idx in range(tag1_idx+1,len(tags))]
        
        combined_result = pd.concat(result, axis=1)
        combined_result["sum"] = combined_result.sum(axis=1)
        combined_result = combined_result[["sum"]]

        resultDataFrame = resultDataFrame.merge(combined_result, left_index = True, right_index = True, how = "outer")
        resultDataFrame.rename(columns = {'sum': ab}, inplace = True)

    assert(resultDataFrame.isnull().values.any() == False)
    #resultDataFrame = resultDataFrame.apply(stats.zscore)
    return(resultDataFrame)


def remove_outliers_from_nbinom(data):

    #remove single AB count values for those that obviously are outside of nbinom distribution
    data_1 = data[(data["ab_id"] == "ITGA6") & (data["ab_count"] < 100)]
    data_2 = data[(data["ab_id"] == "ITGB1") & (data["ab_count"] < 100)]
    data_3 = data[(data["ab_id"] == "TGM1") & (data["ab_count"] < 40)]
    data_4 = data[(data["ab_id"] == "RNApolII") & (data["ab_count"] < 20)]
    data_5 = data[(data["ab_id"] == "Actin") & (data["ab_count"] < 20)]
    data_filter = data_1
    data_filter = data_filter.append(data_2)
    data_filter = data_filter.append(data_3)
    data_filter = data_filter.append(data_4)
    data_filter = data_filter.append(data_5)

    #remove those samples where more than 10 of the 44 ABs are outliers
    occurences_of_samples=data_filter.groupby(['sample_id']).size().to_frame('size')
    occurences_of_samples = occurences_of_samples.reset_index()
    data_filter = data_filter.merge(occurences_of_samples)
    data_filter = data_filter[data_filter["size"] > 10]

    data = data.loc[data['sample_id'].isin(data_filter["sample_id"]) == False]

    return(data)

""" RULES FOR DATA preparation:
    - keep only columns cell==1 (wells with real cells, not empty ones)
    - remove AB with nr==6, bcs it is not measured in half of the data and seems to be unmeasured
    - remove all samples with more than 10 ABs that are obviously out of the nbinom distribution
       (in the histograms there r very obvious AB population with 'zero counts', according to this a threshold was manually selected
        to destinguish an AB count as 'outlier')  
"""
def prepare_data(data):
    data["sample_id"] = data["Barcode_2"] + data["sample_folder"]
    data.rename(columns = {'antibody_count': 'ab_count'}, inplace = True)
    data['sample_id'].astype(str)
    data = data.rename(columns = {"Ab_name" : "ab_id"})

    #remove AB_nr 6 (only measure in around half of all samples)
    data = data[data.Ab_barcode_nr != 6]
    #keep only datapoints that actually come from a non-empty well
    data = data[data["cell"] == 1]
    #remove the samples with more than 10 ABs out of nbinom distribution
    data = remove_outliers_from_nbinom(data)
    
    #calculate relative AB count
    rel_data=data.groupby(['sample_id'])['ab_count'].sum().reset_index()
    rel_data.rename(columns = {'ab_count': 'ab_sample_sum'}, inplace = True)
    rel_data['sample_id'].astype(str)

    data = data.merge(rel_data)
    data["ab_count_relative"] = data.ab_count / data.ab_sample_sum
    data = data[["ab_count", "cell", "ab_id", "sample_id", "ab_count_relative", "Ab_barcode_nr"]]

    return(data)

def main():

    data_file = sys.argv[1]
    data = pd.read_csv(data_file, sep='\t')
    data = prepare_data(data)
    distaneMatrix = distance_from_correlation(data)

    distaneMatrix.to_csv("/Users/t.stohn/Desktop/ab_tag_outlier.tsv")



if __name__ == '__main__':
    main()