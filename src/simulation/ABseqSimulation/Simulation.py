from dataclasses import replace
from turtle import shape
from scipy.stats import nbinom
import re
import numpy as np
from collections import namedtuple
import random
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn': to disable the warning message when subsetting data ('A value is trying to be set on a copy of a slice from a DataFrame.')
import regex as re
from dataclasses import dataclass
import sys
import os
import ast
sys.path.append('./src/methods/ToolBox')
from functions import *
from scipy.linalg import eigh, cholesky
from copulas.multivariate import GaussianMultivariate
from NegativeBinomialCopulaWrapper import NegBinomUnivariate
import copulas
from copulas.univariate.base import BoundedType, ParametricType, ScipyModel
from copulas.multivariate import Multivariate
sys.path.append('./src/methods/nearest_correlation')
from nearest_correlation import nearcorr

def LINE():
    return sys._getframe(1).f_lineno

def replace_negatives(x):
    if isinstance(x, (int, float)) and x < 0:
        return 0
    else:
        return x
    
#data from paper
    # 25 cycles of PCR

    #UMI DISTRIBUTION
    #size 3.83 && mu 4.68

    #PROTEIN DISTRIBUTION
    #most Abcounts are between 100 and 500, with a few beeing greater, and IGTB1 beeing around 23000 // 
    # size parameter has a median of 3.5 beeing mostly around 2.5 and 4

#introduce a biological variance matrix: e.g. neg. binomial distributed
#introduce a technical variance matrix: AB amount differences, annealing differences

def convert_neg_binom_params(mu, size):
    p = size/(size+mu)
    n = mu*p/(1.0 - p)
    return n, p
    
class ProteinCountDistribution():

    def __convert_params(self, mu, size):
        p = size/(size+mu)
        n = mu*p/(1.0 - p)
        return n, p

    def __init__(self, number, mu, size):
        #calculate n and p parameters
        n, p = self.__convert_params(mu, size)
        abCountVector = nbinom.rvs(n ,p, size = number)
        self.abCountVector = abCountVector

    def distributionValues(self):
        return(self.abCountVector)

#generates a ground truth matrix for protein counts
class Parameters():

    #incubation time is considered constant for now
    #antibody count is not considered for now: in sc IDseq constant AB concentration was used for all targets (0.1myg/ml over night)
    #biological variation matrix
    #technical variation matrix
    
    #data types to hold proteinlevels/dist/correlations which are inserted in dictionaries mapping them to individual proteins later on
    ProteinMeanDist = namedtuple('ProteinLevel', ['mu', 'size', 'number'])
    ProteinDist = namedtuple('proteinDist',['mu','size'])
    ProteinCorrelation = namedtuple('ProteinCorrelation', ['mean', 'std', 'number'])

    """ SIMULATION Parameters """
    simulationName = ""
    abBindingEfficiency = None
    seqAmplificationEfficiency = None
    #list of namedtuples ProteinLevel
    ProteinLevels = None
    size = None
    CellNumber = None
    ProteinNumber = None
    abDuplicates = 1
    abDuplicateRange = [1,1]
    abDuplicateDisturbance=0
    pcrCycles=1
    pcrCapture=1
    libSize=[1,1]
    batchFactors=None
    noise=0.0
    noiseExtrinsic=0.0
    proteinNoise=0.0
    
    #CORRELATIONS
    proteinCorrelationDists = []
    proteinCorrelationMean = 0.0
    proteinCorrelationStd = 0.4
    numberProteinCorrelations = 0

    #CLUSTER
    numberOfClusters = 1
    abundanceFactors=[]
    numberClusterSpecificProteins=[]

    correlationFactors=[]
    correlationSets=[]

    cellPercentages=[]

    cellABCountAtGroundtruth = {}

    """ PARAMETERS
        rangeVector: a vector of quadruples(range, number, mean, std) for Abs of different distributions
        abDuplicates: to simulate the usage of several tags for one AB
    """

    def __parseParameters(self, paramFile):
        self.simulationName = os.path.basename(removesuffix(paramFile, '.ini'))
        file = open(paramFile, "r")
        line = file.readline()
        #count the number of proteins/ correlated proteins, they must match in the end
        proteinNumber = 0
        corrNumber = 0
        while line:

            if( not(line.startswith("#")) and ("INIRANGE" in line) ):
                print("ERROR: INI-FILE contains an INIRANGE-variable, those variables can only be used for BENCHMARK INI-FILES to create a range of values, that are then written in its own INI file \
                      for the simulation. Here, we are in the simulation functionality, which only takes INI-FILES with fixed parameters (no ranges!!)")
                exit()
            if(str.startswith(line, "ProteinLevels")):
                info = re.match(("ProteinLevels=\[(.*)\]"), line)
                info = str(info[1]).split(";")

                for element in info:
                    parts = re.match(("\s*\[\s*\[\s*([0-9]*)\s*,\s*(\d+[\.\d+]*)\s*\]\s*,\s*([0-9]*)\s*\]\s*"), element)
                    newProteinLevels = self.ProteinMeanDist(int(parts[1]), float(parts[2]), int(parts[3]))
                    proteinNumber += newProteinLevels.number
                    if(self.ProteinLevels is not None):
                        self.ProteinLevels.append(newProteinLevels)
                    else:
                        self.ProteinLevels = [newProteinLevels]
                self.ProteinNumber = proteinNumber
            elif(str.startswith(line, "size=")):
                info = re.match(("size=(.*)"), line)
                self.size = info[1].rstrip("\n")
            elif(str.startswith(line, "CellNumber=")):
                info = re.match(("CellNumber=(.*)"), line)
                self.CellNumber = int(info[1].rstrip("\n"))
            elif(str.startswith(line, "abDuplicates=")):
                info = re.match(("abDuplicates=(.*)"), line)
                self.abDuplicates = int(info[1].rstrip("\n"))
            elif(str.startswith(line, "abBindingEfficiency=")):
                info = re.match(("abBindingEfficiency=(.*)"), line)
                self.abBindingEfficiency = float(info[1].rstrip("\n"))
            elif(str.startswith(line, "seqAmplificationEfficiency=")):
                info = re.match(("seqAmplificationEfficiency=(.*)"), line)
                self.seqAmplificationEfficiency = float(info[1].rstrip("\n"))
            elif(str.startswith(line, "pcrCycles=")):
                info = re.match(("pcrCycles=(.*)"), line)
                self.pcrCycles = float(info[1].rstrip("\n"))
            elif(str.startswith(line, "pcrCapture=")):
                info = re.match(("pcrCapture=(.*)"), line)
                self.pcrCapture = float(info[1].rstrip("\n"))  
            elif(str.startswith(line, "batchFactors=")):
                info = re.match(("batchFactors=(.*)"), line)
                info = str(info[1]).split(",")
                for batchNum in info:
                    num = float(batchNum)
                    if(self.batchFactors is not None):
                        self.batchFactors.append(num)
                    else:
                        self.batchFactors = [num]   
            elif(str.startswith(line, "libSize=")):
                info = re.match(("libSize=\[(.*)\]"), line)
                info = str(info[1]).split(",")
                assert(len(info)==2)
                self.libSize[0]=float(info[0])
                self.libSize[1]=float(info[1])
            elif(str.startswith(line, "abDuplicateRange=")):
                info = re.match(("abDuplicateRange=\[(.*)\]"), line)
                info = str(info[1]).split(",")
                assert(len(info)==2)
                self.abDuplicateRange[0]=float(info[0])
                self.abDuplicateRange[1]=float(info[1])
            elif(str.startswith(line, "abDuplicateDisturbance=")):
                info = re.match(("abDuplicateDisturbance=(.*)"), line)
                self.abDuplicateDisturbance=float(info[1])
            elif(str.startswith(line, "noise=")):
                info = re.match(("noise=(.*)"), line)
                self.noise = float(info[1].rstrip("\n"))
            elif(str.startswith(line, "proteinNoise=")):
                info = re.match(("proteinNoise=(.*)"), line)
                self.proteinNoise = float(info[1].rstrip("\n"))
            #CORRELATION VARIABLES
            elif(str.startswith(line, "proteinCorrelationDists=")): 
                #same for as ProteinLevels: [[mu, sd, #proteins]; ...] 
                # with as many Correlation Distributions as wanted
                info = re.match(("proteinCorrelationDists=\[(.*)\]"), line)
                info = str(info[1]).split(";")
                for element in info:
                    parts = re.match(("\s*\[\s*(-*\d+[\.\d+]*)\s*,\s*(\d+[\.\d+]*)\s*,\s*([0-9]*)\s*\]\s*"), element)
                    newCorrDist = self.ProteinCorrelation(float(parts[1]), float(parts[2]), int(parts[3]))
                    corrNumber += newCorrDist.number
                    if(self.proteinCorrelationDists is not None):
                        self.proteinCorrelationDists.append(newCorrDist)
                    else:
                        self.proteinCorrelationDists = [newCorrDist]
                assert corrNumber <= self.ProteinNumber, "The numbers of correlated proteins in <proteinCorrelations> can not be higher than proteins in <ProteinLevels>h!! It can be smaller, those proteins will be correlated by drawing form a gaussian around zero (see params)"
            elif(str.startswith(line, "proteinCorrelationMean=")):
                info = re.match(("proteinCorrelationMean=(.*)"), line)
                self.proteinCorrelationMean = float(info[1].rstrip("\n"))
            elif(str.startswith(line, "proteinCorrelationStd=")):
                info = re.match(("proteinCorrelationStd=(.*)"), line)
                self.proteinCorrelationStd = float(info[1].rstrip("\n"))
            #CLUSTER VARIABLES
            elif(str.startswith(line, "abundanceFactors=")):
                info = re.match(("abundanceFactors=(.*)"), line.rstrip('\n'))
                info = info[1].rstrip("\n")
                assert info[0] == "[", "First character must be '[' for following line in the ini file: " + line
                assert info[-1] == "]", "Last character must be ']' for following line in the ini file: " + line
                self.abundanceFactors = ast.literal_eval(info)
                self.abundanceFactors = [float(element) for element in self.abundanceFactors]
            elif(str.startswith(line, "numberClusterSpecificProteins=")):
                info = re.match(("numberClusterSpecificProteins=(.*)"), line.rstrip('\n'))
                info = info[1].rstrip("\n")
                assert info[0] == "[", "First character must be '[' for following line in the ini file: " + line
                assert info[-1] == "]", "Last character must be ']' for following line in the ini file: " + line
                self.numberClusterSpecificProteins = ast.literal_eval(info)
                self.numberClusterSpecificProteins = [int(element) for element in self.numberClusterSpecificProteins]
            elif(str.startswith(line, "correlationFactors=")):
                info = re.match(("correlationFactors=(.*)"), line.rstrip('\n'))
                info = info[1].rstrip("\n")
                assert info[0] == "[", "First character must be '[' for following line in the ini file: " + line
                assert info[-1] == "]", "Last character must be ']' for following line in the ini file: " + line                
                self.correlationFactors = ast.literal_eval(info)
                self.correlationFactors = [float(element) for element in self.correlationFactors]
            elif(str.startswith(line, "cellPercentages=")):
                info = re.match(("cellPercentages=(.*)"), line.rstrip('\n'))
                info = info[1].rstrip("\n")
                assert info[0] == "[", "First character must be '[' for following line in the ini file: " + line
                assert info[-1] == "]", "Last character must be ']' for following line in the ini file: " + line
                self.cellPercentages = ast.literal_eval(info)
                self.cellPercentages = [float(element) for element in self.cellPercentages]
            elif(str.startswith(line, "correlationSets=")): 
                #same for as ProteinLevels: [[mu, sd, #proteins]; ...] 
                # with as many Correlation Distributions as wanted
                info = re.match(("correlationSets=\[(.*)\]"), line)
                info = str(info[1]).split(";")
                for element in info:
                    newSet = ast.literal_eval(element)
                    newSet = [int(element) for element in newSet]
                    self.correlationSets.append(newSet)
            line = file.readline()
        #we need exactly one abundance/ correlation/ percentage factor for every clusters
        self.numberOfClusters = len(self.cellPercentages)
        #set cellPercentages to 100 if no clusters exist
        if(self.cellPercentages == []):
            self.cellPercentages = [100]
            self.numberOfClusters = 1
            
        assert ( len(self.abundanceFactors) + 1   == self.numberOfClusters), "Error in ini file: cellPercentages and abundanceFactors are of different length (we need ONE FACTOR per ADDITIONAL CLUSTER)"
        assert ( len(self.numberClusterSpecificProteins) + 1   == self.numberOfClusters), "Error in ini file: cellPercentages and number of scaled proteins per cluster <numberProteins> are of different length (we need ONE FACTOR per ADDITIONAL CLUSTER)"
        
        #we also have cluster specific correlation variables
        if( len(self.correlationFactors) > 0 or len(self.correlationSets) > 0):
            assert ( len(self.correlationFactors) + 1 == self.numberOfClusters), "Error in ini file: cellPercentages and correlationFactors are of different length (we need ONE FACTOR per ADDITIONAL CLUSTER)"
            assert ( len(self.correlationSets) + 1 == self.numberOfClusters), "Error in ini file: cellPercentages and correlationSets are of different length (we need ONE set of diff-corr protein distributions per ADDITIONAL CLUSTER)"

        assert (sum(self.cellPercentages) == 100), "Sum of cellPercentages must be 100"
        if(self.correlationSets is not None):
            for CorrSet in self.correlationSets:
                assert len(CorrSet)<=len(self.proteinCorrelationDists), "CorrelationSet Error: per cluster we can only scale as many correlationSets as there are distirbutions of correlations to sample from..."
        totalNumberIntercorrelatedProteins = 0
        for correlationDist in self.proteinCorrelationDists:
            totalNumberIntercorrelatedProteins += correlationDist.number
        assert(totalNumberIntercorrelatedProteins <= self.ProteinNumber), "The number of protein-protein correlation stated in the 3rd parameters in proteinCorrelationDists does not sum up to the number of proteins measured, we have more correlations than possible!!!"
            
    def __init__(self, paramter_file):
        self.__parseParameters(paramter_file)

class SingleCellSimulation():

    """ PARAMETERS """
    parameters = None
    ab_sampled = None
    output_dir = "./bin/SIMULATIONS/"
    proteinDistributions = {} #the actual distributions for every protein
    #(this is different to proteinDist from parameters, which stores the distributions of mu, size for the whole dataset)

    """ MODELLED DATA """
    groundTruthData = None
    #sample * AB matrix of simulated data
    simulatedData = None

    def __init__(self, parameters):
        self.parameters = parameters

    """ we see a strong dependance of size on mu for single-cell protein distributions, therefore calualte size depending on mu"""
    def __calculate_size_from_mu(self, mu):
        size_formula_with_mu_insertion = self.parameters.size.replace('MU', str(mu))
        result = eval(size_formula_with_mu_insertion) 
        return(result)
        
    """ functions to simulate the cell * proteinCount matrix"""
    def __generate_neg_binom_distributed_protein_counts(self):
        proteinCount = 1
        proteinCountMatrix = None
        for proteinRange in self.parameters.ProteinLevels:
            for i in range(proteinRange.number):
                #for every protein simulate a neg.binom distibuted number for every cell
                mu = random.randrange(proteinRange.start, proteinRange.end)
                #TO DO: size variance of neg.binom distribution is simpy a random value +- 1
                size_variance = random.randrange(-1, 1)
                size = self.parameters.size + size_variance
                dist = ProteinCountDistribution(self.parameters.CellNumber, mu, size)
                proteinCountVector = dist.distributionValues()
                abName = "AB" + str(proteinCount)
                #add it to a matrix as a new column
                if(proteinCountMatrix is None):
                    proteinCountMatrix = pd.DataFrame({abName : proteinCountVector}) 
                    proteinCountMatrix.index = ["sample_" + str(j+1) for j in range(self.parameters.CellNumber)]
                else:
                    proteinCountMatrix[abName] = proteinCountVector

                proteinCount +=1

                proteinDist = self.parameters.ProteinDist(int(mu), float(size))
                self.proteinDistributions[abName] = proteinDist
        
        self.cellABCountAtGroundtruth = proteinCountMatrix.sum(axis = 1)

        proteinCountMatrix = proteinCountMatrix.reset_index().rename(columns={ 'index' : 'sample_id'})
        proteinCountMatrix = proteinCountMatrix.melt(id_vars = ["sample_id"], var_name="ab_id", value_name="ab_count")
        return(proteinCountMatrix)

    """ model treatment effect into data """
    def __insert_treatment_effect(self, data):

        result = []
        concatedTreatments = None
        
        #generate lists of sampleIds for each treatment
        sampleIdVector = (data["sample_id"].unique())
        sampleIdList = np.array_split(sampleIdVector, len(self.parameters.treatmentVector)+1 )
        print("Simulating " + str(len(self.parameters.treatmentVector)) + " treaments.")
        i = 0
        
        #model each treatment with its own samples (except for the first one, that one stays uneffected)
        dataTreatment = data[data["sample_id"].isin(sampleIdList[0])]
        #IMPORTANT: control treatment is called control, this is hard coded
        #nad used again for benchmarking
        dataTreatment["cluster_id"] = "control"
        result.append(dataTreatment)
        for sampleIds in sampleIdList[1:]:
            #SHIFTVALUES IS VECTOR OF FACTORS AND DATATREATMENT IS VECTOR OF SAMPLE IDS
            dataTreatment = data[data["sample_id"].isin(sampleIds)]
            shiftValues = self.parameters.treatmentVector[i]
            #for each differentially expressed protein
            j = 0
            for value in shiftValues:
                proteinId = self.parameters.diffExProteins[i][j]
                dataTreatment.loc[dataTreatment["ab_id"] == proteinId, "ab_count"] *= value
                j+=1
            #IMPORTANT: treatments are enumerated from 0 to x, this is later again used when
            #benchmarking the treatment effect
            dataTreatment["cluster_id"] = "cluster_" + str(i)
            result.append(dataTreatment)
            i+=1
        concatedTreatments = pd.concat(result)
        self.groundTruthData = concatedTreatments

        return(concatedTreatments)

    def __rescale_AB_counts(self, data):

        #calculate new total AB count per sample
        newToalABCount = {}
        for sampleID in (data["sample_id"]).unique():
            newToalABCount[sampleID] = data.loc[data['sample_id'] == sampleID, 'ab_count'].sum()

        #for every AB of this sample, rescale by factor diff
        for sample in (data["sample_id"]).unique():
            #calculate diff between new and od count
            diff = self.cellABCountAtGroundtruth[sample] / newToalABCount[sample]

            for ab in (data["ab_id"]).unique():
                data.loc[data["sample_id"] == sample, "ab_count"] *= diff

        #round everything down
        data['ab_count'] = data['ab_count'].apply(np.round)
        return(data)

    """ model batch effect """
    def __insert_batch_effect(self, data):
        batches = len(self.parameters.batchFactors)
        result = []

        #generate lists of sampleIds for each batch
        sampleIdVector = (data["sample_id"].unique())
        sampleIdList = np.array_split(sampleIdVector, batches)
        print("Simulating " + str(batches) + " batches.")

        #generate a new dataframe for each batch
        i = 0
        for sampleIds in sampleIdList:
            dataBatch = data[data["sample_id"].isin(sampleIds)]
            dataBatch["ab_count"] *= self.parameters.batchFactors[i]
            dataBatch["batch_id"] = "batch_" + str(i)
            result.append(dataBatch)
            i+=1

        concatedBatches = pd.concat(result)
        return(concatedBatches)

    """ model libsize effect """
    def __insert_libsize_effect(self, data):
        print("Simulating different libsizes.")
        
        perturbFrame = pd.DataFrame()
        perturbFrame["sample_id"] = data["sample_id"].unique()
        #sample lbisize factors from normal dist: 1.value = mean, 2.= std
        perturbFrame["factor"] = np.random.lognormal(self.parameters.libSize[0],self.parameters.libSize[1], len(perturbFrame))

        perturbDict = dict(zip(perturbFrame["sample_id"], perturbFrame["factor"]))
        data["ab_count"] = data["ab_count"] * data["sample_id"].map(perturbDict)

        return(data)

    def __covariance(self, x, y, corr):
        
        sdX = np.std(x)
        sdY = np.std(y)
        sdProduct = sdX * sdY
        
        cov = corr * sdProduct
        return cov

    #generate correlated counts by the help of the cholevski method (outdated version without correlations)
    def __correlate_proteins_cholevski(self, data, prot1, prot2, corr):
        
        #generate multivariate gaussian dist data
        x_previous = np.array(data.loc[data["ab_id"] == prot1,"ab_count"])
        y_previous = np.array(data.loc[data["ab_id"] == prot2,"ab_count"])
        
        x_previous = np.random.normal(np.mean(x_previous), np.std(x_previous), len(x_previous))
        y_previous = np.random.normal(np.mean(y_previous), np.std(y_previous), len(y_previous))
        covM = np.array([
                            [np.var(x_previous), self.__covariance(x_previous, y_previous, corr)],
                            [self.__covariance(x_previous, y_previous, corr), np.var(y_previous)]
                        ], dtype = float)
        means = [x_previous.mean(), y_previous.mean()]  
        values = np.random.multivariate_normal(means, covM, len(x_previous)).T

        #generating negbinom and map to same idx as gaussian data
        #keeoing rank for correlation but skewing data by negbinom
        dataFrame = pd.DataFrame({'X':values[0], 'Y':values[1]})
        dataFrame['Idx'] = np.arange(len(dataFrame))
        
        dist1 = self.proteinDistributions[prot1]
        dist2 = self.proteinDistributions[prot2]
        distX = ProteinCountDistribution(self.parameters.CellNumber, dist1.mu, dist1.size)
        distY = ProteinCountDistribution(self.parameters.CellNumber, dist2.mu, dist2.size)
        proteinCountVectorX = distX.distributionValues()
        proteinCountVectorY = distY.distributionValues()
        negBinomX = np.sort(proteinCountVectorX)
        negBinomY = np.sort(proteinCountVectorY)

        dataFrameSortX = dataFrame.sort_values(by=['X'])
        dataFrameSortX['Xnb'] = negBinomX
        dataFrameSortY = dataFrameSortX.sort_values(by=['Y'])
        dataFrameSortY['Ynb'] = negBinomY

        newData = dataFrameSortY.sort_values(by=['Idx'])
        returnData = np.array([newData['Xnb'], newData['Ynb']])

        return(returnData)
    
    #creates sets of proteins that are CORRELATED: these sets do not contain overlapping proteins
    def __groupProteinsForCorrelationSets(self):
        correlationProteinsSets = []
        proteinList = range(self.parameters.ProteinNumber)
        for proteinCorr in self.parameters.proteinCorrelationDists:
            subset = random.sample(proteinList, proteinCorr.number)
            correlationProteinsSets.append(subset)
            #remove elements that were already selected for next sampling step
            proteinList = [element for element in proteinList if element not in subset]
        return(correlationProteinsSets)
    
    def __sampleFromCorrelationDist(self, setIdx=None, scalingFactor = 1.0):
        mean = 0.0
        std = 0.0
        if(setIdx is None):
            mean = self.parameters.proteinCorrelationMean * scalingFactor
            std = self.parameters.proteinCorrelationStd
            print(self.parameters.proteinCorrelationStd)
            print(str(mean))

        else:
            mean = self.parameters.proteinCorrelationDists[setIdx].mean * scalingFactor
            std = self.parameters.proteinCorrelationDists[setIdx].std

        corr = np.random.normal(mean, std)
        return(corr)

    def __isCorrelationPair(self, correlationProteinsSets, i, j):
        #if the pair is in any of the correlation sets
        for setIdx in range(len(self.parameters.proteinCorrelationDists)):
            #if the pair is in a Set draw from that dist
            if( (i in correlationProteinsSets[setIdx]) and 
                (j in correlationProteinsSets[setIdx])):
                    return setIdx
        return -1            
                            
    def __isScaledPair(self, setIdx, scaledSetsPerCluster, clusterIdx):                      
        if(len(scaledSetsPerCluster) > 1 and setIdx in scaledSetsPerCluster[clusterIdx]):
            assert clusterIdx!=0, "baseline cluster SHOULD NOT GET CORRELATIONS SCALED!!"
            return(True)
        return(False)
    
    #method according to Levandowski et. al 2009
    # better than second method (from post see below) since it also introduces neg. correlations
    # implementation based on matlab code from https://stats.stackexchange.com/questions/124538/how-to-generate-a-large-full-rank-random-correlation-matrix-with-some-strong-cor
    # 0.2 - 10 is a 'good' parameter range for strong to weak correlations for 60 proteins (manually observed)
    def __vineBeta(d, betaparam):
        P = np.zeros((d, d))           # storing partial correlations
        S = np.eye(d)

        for k in range(1, d):
            for i in range(k + 1, d):
                P[k, i] = np.random.beta(betaparam, betaparam) # sampling from beta
                P[k, i] = (P[k, i] - 0.5) * 2     # linearly shifting to [-1, 1]
                p = P[k, i]
                for l in range(k - 1, 0, -1): # converting partial correlation to raw correlation
                    p = p * np.sqrt((1 - P[l, i] ** 2) * (1 - P[l, k] ** 2)) + P[l, i] * P[l, k]
                S[k, i] = p
                S[i, k] = p

        # permuting the variables to make the distribution permutation-invariant
        permutation = np.random.permutation(d)
        S = S[permutation][:, permutation]
        return S

    # 1.) correlationProteinsSets: list of lists, where every sublist contains all proteins that r intercorrelated/ order of proteinCorrelationDist
    # 2.a) scaledSetsPerCluster: list of lists, where every sublist contains which correlationproteinSets should be scaled: ordered by cluster (first element zero for baseline)
    # 2.b) correlationScalingFactors: list of all the scaling factors for the proteinSets that r supposed to be scaled: ordered by clusters (first element zero for baseline)
    # 3.) clusterIdx: current cluster with 0==Baseline
    def __generateCovarianceMatrix(self, correlationProteinsSets, scaledSetsPerCluster, correlationScalingFactors, clusterIdx):

        #initialize zero covariance matrix
        covMat = [[0.0 for _ in range(self.parameters.ProteinNumber)] for _ in range(self.parameters.ProteinNumber)]
        
        #iterate through all protein-pairs (triangular matrix - but fill all values)
        for i in range(self.parameters.ProteinNumber):
            isWeighted = False
            for j in range(i, self.parameters.ProteinNumber):
                if(i == j):
                    covMat[i][j] = 1.0
                    covMat[j][i] = 1.0
                else:
                    #check if pair is part of a correlated set
                    #if true returns setIdx, otherwise minus one
                    setIdx = self.__isCorrelationPair(correlationProteinsSets, i, j)
                    isScaled = False
                    if(setIdx >= 0):
                        isScaled = self.__isScaledPair(setIdx, scaledSetsPerCluster, clusterIdx)                  

                    if( (setIdx >= 0) and isScaled):
                        scalingFactor = correlationScalingFactors[clusterIdx]
                        corr = self.__sampleFromCorrelationDist(setIdx, scalingFactor)
                    elif(setIdx >= 0):
                        corr = self.__sampleFromCorrelationDist(setIdx)
                    else:
                        corr = self.__sampleFromCorrelationDist(None)
 
                    if(corr > 1.0): 
                        corr = 1.0
                    elif(corr < -1.0): 
                        corr = -1.0
                    covMat[i][j] = corr
                    covMat[j][i] = corr
                    
        #check if matrix is positive semi definite
        print("___BEORE___\n")
        print(np.round(covMat,2))
        minEigenVal = min(np.linalg.eigvals(covMat))
        if(minEigenVal <0):
            covMat = nearcorr(np.array(covMat), max_iterations=10000)
        print("___AFTER___\n")
        print(np.round(covMat,2))

        return(covMat)
        
    #make a list of all the protein distributions
    def __generateProteinDistributions(self,scaledProteinIdxs, clusterIdx):
            
        distributions = []
        protein = 0
        #we need to insert 1.0 abudnace factor for baseline clusterIdx=0
        abundanceFactorsWithBaseline = self.parameters.abundanceFactors
        abundanceFactorsWithBaseline.insert(0, 1.0) #base cluster has no sclaing of protein abudnance (so times factor 1.0)

        proteinsToScale = scaledProteinIdxs[clusterIdx]
        for proteinRange in self.parameters.ProteinLevels:
            
            #generate mu/ size for all proteins of the proteinRange distribution (sampled from distribution of means and the calculating size)
            dist = ProteinCountDistribution(proteinRange.number, proteinRange.mu, proteinRange.size)
            proteinCountVector = dist.distributionValues() #vector of mean values for all proteins for this proteinRange-distribution
             
            for i in range(proteinRange.number):       
                
                #get the very specific mu/ size for the protein i in this proteinRange-distribution
                mu_tmp = proteinCountVector[i]
                size_tmp = self.__calculate_size_from_mu(mu_tmp)
                         
                #sclae mu by cluster specific scaling factor
                if(protein in proteinsToScale):
                    assert clusterIdx != 0, "The zero Cluster(baseline, not sclaed) can not have a scaling factor"
                    mu_tmp = mu_tmp * abundanceFactorsWithBaseline[clusterIdx]

                n,p = convert_neg_binom_params(mu_tmp, size_tmp)
                tmpDict = {'loc': 0.0, 'n': n, 'p': p, 'type': NegBinomUnivariate}
                distributions.append(tmpDict)

                protein += 1

        return(distributions)
    
    def __sampleCellTimesProteinMatrix(self, percentage, dist, cov, clusterIdx):

        print(np.round(cov,2))
        print(dist)
        
        copula_params = {}
        copula_params['correlation'] = cov
        copula_params['univariates'] = dist
        copula_params['type'] = copulas.multivariate.gaussian.GaussianMultivariate
        copula_params['columns'] = ["AB" + str(number) for number in range(self.parameters.ProteinNumber)]        
        copula_model = Multivariate.from_dict(copula_params)

        cellNumber = round(self.parameters.CellNumber * (percentage/100.0))
        copulaResult = copula_model.sample(cellNumber)

        #bring copulaResult in the right format
        copulaResult = copulaResult.reset_index().rename(columns={ 'index' : 'sample_id'})
        copulaResult = copulaResult.melt(id_vars = ["sample_id"], var_name="ab_id", value_name="ab_count")
        copulaResult.insert(0, 'batch_id', 'batch')
        copulaResult.insert(0, 'ab_type', 'any_ab')

        if(clusterIdx == 0):
            copulaResult.insert(0, 'cluster_id', 'control')
        else:
            copulaResult.insert(0, 'cluster_id', ('cluster_' + str(clusterIdx)))

        return(copulaResult)
            
         
    def __generateClusterSpecificCovarianceMatrices(self):
    
        covariancematrix = []
        
        # list of protein Idx for every dist in proteinCorrelationDists
        # we have the number of intercorrelated proteins from proteinCorrelationDist variable
        # now we sample this number fo proteins for every proteinDist that is stated in the variable
        correlationProteinsSets = []
        #stores a list of lists: lists for every proteinCorrelationDist with a list of all proteins belonging to this dist
        correlationProteinsSets = self.__groupProteinsForCorrelationSets()
        
        print(correlationProteinsSets)
        #proteinSet from proteinCorrelationDist for which CORRELATIONS are SCALED differently
        scaledCorrelationIdxs = self.parameters.correlationSets
        scaledCorrelationIdxs.insert(0, []) #base cluster has no sclaing of correlations
        
        #scaling factors for all correlations in the scaledCorrelationIdxs
        correlationScalingFactors = self.parameters.correlationFactors
        correlationScalingFactors.insert(0, 1.0)
        for idx in range(self.parameters.numberOfClusters):
            covariancematrix.append(self.__generateCovarianceMatrix(
                correlationProteinsSets, scaledCorrelationIdxs, correlationScalingFactors, idx))

        return(covariancematrix)
            
    def __generateGroundTruthWithCopula(self):
        # ALL LIST ORDERED BY CLUSTER
        distributions = [] #list of the cluster specific distribution for proteins (scaled by factors for clusters)
        #it is a list of lists: for every cluster, for every proteins, the neg.binom. parametrization
        groundTruthDataList = []

        #SAMPLE PROTEIN IDs that are cluster specific (add empty set for baseclass = zero scaled proteins)
        clusterSpecificProteinIdxs = [[]]
        if(self.parameters.numberOfClusters > 1):
            #idices go only from 0 to clsuters-1, since we do not need to get clusterSpecificProtein nuber for baseline cluster,
            #and the numberClusterSpecificProteins list is ordered starting from 0 with the first cluster, so that zero does not refer to baseline here
            for idx in range(self.parameters.numberOfClusters-1):
                subset = random.sample(range(self.parameters.ProteinNumber), self.parameters.numberClusterSpecificProteins[idx])
                clusterSpecificProteinIdxs.append(subset)
                
        #CLUSTER SPECIFIC PROTEIN-DISTRIBUTION-PARAMETRIZATION (scaled for proteins)
        for idx in range(self.parameters.numberOfClusters):
            distributions.append(self.__generateProteinDistributions(clusterSpecificProteinIdxs, idx))

        #CLUSTER SPECIFIC COVARIANCE MATRIX (scaled for protein-apirs) - list: one per cluster
        covariancematrix =  self.__generateClusterSpecificCovarianceMatrices()
        
        #sample data from copula for every cluster and combine
        for idx in range(self.parameters.numberOfClusters):
            groundTruthDataList.append(self.__sampleCellTimesProteinMatrix(self.parameters.cellPercentages[idx], 
                                                distributions[idx],
                                                covariancematrix[idx],
                                                idx))

        combined_df = pd.concat(groundTruthDataList, ignore_index=True)
        return(combined_df)

    """ model correlated proteins """
    #this is without introducing correlations (outdated version)
    def __insert_correlations_between_proteins_cholevski(self, data):
        #sort data first according to sample
        #(when we change values for certain proteins we add values of another protein, the returned vector
        # might not be in the same order of sampels if not sorted first)
        data = data.sort_values(by=['sample_id'])

        for corr in self.parameters.proteinCorrelations:
            #origional negbinom distribution + factor times dependant protein count
            newCorrelatedValues = self.__correlate_proteins_cholevski(data,
                                                            corr.prot1,
                                                            corr.prot2, 
                                                            corr.factor)
            data.loc[data["ab_id"] == corr.prot1,"ab_count"] = newCorrelatedValues[0,:].copy()
            data.loc[data["ab_id"] == corr.prot2,"ab_count"] = newCorrelatedValues[1,:].copy()

        data.loc[data["ab_count"] < 0,"ab_count"] = 0

        self.groundTruthData = self.groundTruthData.sort_values(by=['sample_id'])
        self.groundTruthData = data
        return(data)

    """ model duplicates for ABs """
    #call this function as last for generating ground truth, otherwise the duplicates could be distributed in different bins for
    #treatment/ protein correlation effects
    def __insert_ab_duplicates(self, data):
        print("Simulating AB duplicates.")
        firstAbData = data.copy()
        firstAbData["ab_id"] = firstAbData["ab_id"] + "a"
        result = [firstAbData]
        abIdVector = (data["ab_id"].unique())
        for abId in abIdVector:
            for i in range(1,self.parameters.abDuplicates):
                n = random.uniform(self.parameters.abDuplicateRange[0],self.parameters.abDuplicateRange[1])
                dataSample = data[data["ab_id"] == abId]
                dataSample["ab_id"] = dataSample["ab_id"] + chr(97+i)
                #Ab duplicate is n-times the origional value
                dataSample["ab_count"] = (n * dataSample["ab_count"])
                #and is disturbed by a normal dist, the variance is variance parameter times average of origional AB count
                dataSample["ab_count"] = np.random.normal(dataSample["ab_count"], self.parameters.abDuplicateDisturbance*np.mean(data[data.ab_id==abId]["ab_count"]))
                dataSample["ab_count"] = dataSample["ab_count"].round(decimals = 0)
                dataSample.loc[dataSample["ab_count"] < 0,"ab_count"] = 0
                result.append(dataSample)
        concatedSamples = pd.concat(result)
        return(concatedSamples)

    """ Generating GroundTruth of the Protein abundancies in all single cells """
    #outdated method without correlations to sample from neg. binom distribution for simulations
    def __generateSimpleGroundTruth(self, parameters):
        #generate matrix of cells * proteinCounts
        proteinCountMatrix = self.__generate_neg_binom_distributed_protein_counts()
        
        proteinCountMatrix.insert(0, 'cluster_id', 'cluster')
        proteinCountMatrix.insert(0, 'batch_id', 'batch')
        proteinCountMatrix.insert(0, 'ab_type', 'any_ab')
        
        self.groundTruthData = proteinCountMatrix

    """ Simulating the Detection of the GroudnTruth Protein Counts """
    def __simulate_ab_binding(self, data):
        nonZeroEntries = np.count_nonzero(data.ab_count)
        number = int(self.parameters.abBindingEfficiency * nonZeroEntries)

        return(data)

    def __generateUmiData(self, data):

        #reset index before repeating index elements
        data = data.reset_index(drop=True)
        umiData = data.loc[data.index.repeat(data["ab_count"])]
        umiData["umi_id"] = range(len(umiData.index))
        return(umiData)

    def __pcrAmplify(self, data, readNum = 0):
    
        pcrNumber = int(self.parameters.pcrCapture * len(data.index))

        for i in range(int(self.parameters.pcrCycles)):
            print("PCR Cycle " + str(i))
            tmp_readsToPcrAmplify = data.sample(n=pcrNumber, replace = False, random_state=1)
            data = pd.concat([data, tmp_readsToPcrAmplify])
        
        #sample from all UMIs and remove umis that occur several times
        print("Sampling from UMIs")
        seqNumber = 0
        if(readNum == 0):
            seqNumber = int(self.parameters.seqAmplificationEfficiency * len(data.index))
        else:
            seqNumber = int(readNum)
        data = data.sample(n=seqNumber, replace = False, random_state=1)
        data = data.drop_duplicates()
        
        #concatenate all UMI reads to one value
        data = data.drop(columns=["umi_id"])
        #concatenate again reads for same Protein in same cell after removing UMI column
        print("Concatenate all reads again to protein counts per cell")
        #count the grouped rows, col after groupby doesn t matter,can be anything (sample_id used here)
        data["ab_count"] = data.groupby(['sample_id', 'ab_id'])['sample_id'].transform('size')

        concatenatedData = data.drop_duplicates()

        #insert zeroes for all the no sampled values 
        # (therefore pivot into wider shape to easily detect nans, 
        # then convert them to zeroes and convert back to origional format)
        #save all the columns that we have to keep in this process (ab_id and ab_count are the ones we use for pivoting)
        idColumns = concatenatedData.columns
        idColumns = list(idColumns.drop(["ab_id","ab_count"]))

        dataPivoted=concatenatedData.pivot(index=idColumns, columns='ab_id')['ab_count'].reset_index()
        dataPivoted.columns.name = None
        dataPivoted = dataPivoted.fillna(0)
        data = dataPivoted.melt(id_vars = idColumns, var_name="ab_id", value_name="ab_count")

        return(data)
        
    def __simulate_sequencing_binding(self, data):
        #simulate UMIs
        #for every line simulate UMIs according to ab_count column
        print("Generate single reads per protein count with UMI. Number data points[cell x AB]: " + str(data.shape[0]))
        umiData = self.__generateUmiData(data)
        print("Simulating with a library size of: " + str(len(umiData)))

        #PCR amplify those reads
        #sampling several times(for each PCR cycle) wo replacement and combine the old and sampled dataset
        # (previously simply sample w replacement and a sampling number > 1: however this neglects the idea that reads sampled in the first PCR
        # cycles are more likely to be overrepresented in the downrun)
        print("Generate PCR amplifications of reads")
        pcrData = self.__pcrAmplify(umiData)

        return(pcrData)

    #SAME AS PREVIOUS VERSION, BUT CELLWISE PCR AMPLIFICATION AND SUBSEQUENT SAMPLING (LOADING READS ON SEQUENCER)
    #on top of that SUBSAMPLING after PCR is not performed on perc but an absolute threshold
    #this leads to no skewing in between cells, in case one cell has a much higher protein than others, 
    #and also leads to possibly new distortions bcs of absolute counts beeing selected, the perc selection
    #might not change the underlying distribution a lot
    #=> in total the two methods seemed to not change a lot, looking at totalABcount and proteinWiseABcount distributions
    def __simulate_sequencing_binding_2(self, data):

        result = []
        #generate lists of sampleIds for each batch
        sampleIdVector = (data["sample_id"].unique())

        #generate a new dataframe for each batch
        i = 0
        readNum = (data['ab_count'].sum() / len(sampleIdVector)) * self.parameters.seqAmplificationEfficiency
        for sampleId in sampleIdVector:

            dataSingleCell = data[data.sample_id == sampleId]
            #for every line simulate UMIs according to ab_count column
            umiData = self.__generateUmiData(dataSingleCell)
            pcrData = self.__pcrAmplify(umiData, readNum)

            result.append(pcrData)

        concatedBatches = pd.concat(result)
        return(concatedBatches)

    def __perturb(self,data):

        printToTerminalOnce("\tAdding random noise to data")

        #pertubation per protein
        perturbFrame = pd.DataFrame()
        perturbFrame["ab_id"] = data["ab_id"].unique()
        perturbFrame["factor"] = np.random.normal(1,self.parameters.proteinNoise,len(perturbFrame))
        perturbDict = dict(zip(perturbFrame["ab_id"], perturbFrame["factor"]))
        data["ab_count"] = data["ab_count"] * data["ab_id"].map(perturbDict)

        #cell instrinsic pertubation
        randomVector = np.random.normal(1,self.parameters.noise,len(data))
        data["ab_count"] *= randomVector

        #between cell pertubation
        #disable for now, since it only correlates the samples
            #perturbFrame = pd.DataFrame()
            #perturbFrame["sample_id"] = data["sample_id"].unique()
            #perturbFrame["factor"] = np.random.normal(1,self.parameters.noiseExtrinsic,len(perturbFrame))
            #perturbDict = dict(zip(perturbFrame["sample_id"], perturbFrame["factor"]))
            #data["ab_count"] = data["ab_count"] * data["sample_id"].map(perturbDict)
        data["ab_count"] = data["ab_count"].round(decimals = 0)

        return(data)

    def __add_batch_effect_to_ground_truth(self, newData):
        self.groundTruthData = pd.merge(self.groundTruthData,newData[['sample_id','ab_id','batch_id']],on=['sample_id','ab_id'], how='left')
        
    """ MAIN FUNCTION: 
    1. generates ground truth & 
    2. simulates the protein count detection
    """
    def simulateData(self):
                
        #for every cluster: scale protein distributions; scale correlations; sample using gaussian copula
        # (protein counts scaled/ correlations scaled/ and then sampled from theirs dists by copula)
        self.groundTruthData = self.__generateGroundTruthWithCopula()

        #not interesting at this point (initialy thought of AB duplicates for normalization)
        if(self.parameters.abDuplicates > 1):
            self.groundTruthData = self.__insert_ab_duplicates(self.groundTruthData)

        #scale counts for batch effects (NOT USED AT THIS POINT)
        perturbedData = self.groundTruthData.copy(deep=True)
        if(not (self.parameters.batchFactors is None)):
            perturbedData = self.__insert_batch_effect(perturbedData)
            self.__add_batch_effect_to_ground_truth(perturbedData)
            
        #scale counts for library size effects
        if(not (self.parameters.libSize[0]==1 and self.parameters.libSize[1]==1)):
            perturbedData = self.__insert_libsize_effect(perturbedData)

        #simulate AB binding efficiency
        #discard a fraction of proteinCounts as no AB has bound to them
        #tmp_simulatedData = self.__simulate_ab_binding(perturbedData)

        #simulate PCR amplification and sequencing
        #sampling with replacement to simulate PCR amplification as well as missing out on reads during washing/ sequencing
        #tmp_simulatedData = self.__simulate_sequencing_binding_2(perturbedData)
        
        #add random noise to simulations
        tmp_simulatedData = self.__perturb(perturbedData)

        #remove negative counts from data, and substitute with zero
        self.simulatedData = tmp_simulatedData
        self.simulatedData = self.simulatedData.applymap(replace_negatives)

        return(self.simulatedData)

    """ SAVE THE DATA """
    def save_data(self):
        #safe data
        printToTerminalOnce("\tSave Data\n")

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

        
        self.simulatedData.to_csv(self.output_dir + "/" + self.parameters.simulationName + "_SIMULATED.tsv", sep='\t', index = False)
        self.groundTruthData.to_csv(self.output_dir + "/" + self.parameters.simulationName + "_GROUNDTRUTH.tsv", sep='\t', index = False)
        printToTerminalOnce("\tData saved\n")

