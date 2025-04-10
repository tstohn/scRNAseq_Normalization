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
import copy
import math

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
    #n = number of successes/ p = probability
    #p = size/(size+mu) - same as formular below
    # formular for probability here: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/NegBinomial.html    
    #in nbinom (python) n is number of successes, in R this is already size, n is the number of observations
    
    n = size
    p = 1 / (1 + (mu / size))

    #n = (mu * p)/(1 - p)

    return n, p
    
#class to handle the mu s for all proteins
#these are smapled from another distribution/ and have a Min cutoff
class ProteinCountDistribution():

    def __convert_params(self, mu, size):
        #p = size/(size+mu)
        
        n = size
        p = 1 / (1 + (mu / size))
    
        #n = (mu * p)/(1 - p)

        return n, p

    def __init__(self, number, mu, size, MinProtein=0):
        
        #calculate n and p parameters
        n, p = self.__convert_params(mu, size)
        #instead of drawing number samples, draw one by one and check if theyr above threshold
        #abCountVector = nbinom.rvs(n ,p, size = number)
        abCountVector=np.array([], dtype=int)
        while len(abCountVector) < number:
            draw = nbinom.rvs(n=n, p=p, size=1)
            if draw >= MinProtein:
                abCountVector = np.concatenate((abCountVector, draw))
        self.abCountVector = abCountVector
        
        print("MADE PROTEIN MEANS:")
        print(self.abCountVector)
        
    def distributionValues(self):
        return(self.abCountVector)

#generates a ground truth matrix for protein counts
class Parameters():

    #incubation time is considered constant for now
    #antibody count is not considered for now: in sc IDseq constant AB concentration was used for all targets (0.1myg/ml over night)
    #biological variation matrix
    #technical variation matrix
    
    #store its own file name (needed to write benchmarkIniFile)
    filePath = ""
    
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
    MinProtein = 0
    size = None
    CellNumber = None
    ProteinNumber = None
    abDuplicates = 1
    abDuplicateRange = [1,1]
    abDuplicateDisturbance=0
    pcrCycles=1
    pcrCapture=1
    libSize=[1.0,0.0] #default lib-size factor is 1.0 introducing NO lib-size effect
    batchFactors=None
    noise=0.0
    noiseExtrinsic=0.0
    proteinNoise=0.0
    
    #CORRELATIONS 1.)
    randomProteinCorrelations=0
    betaParameter=0
    addBetaParameters=[]
    
    #CORRELATIONS 2.)
    proteinCorrelationDists = []
    proteinCorrelationMean = 0.0
    proteinCorrelationStd = 0.4
    numberProteinCorrelations = 0

    #CLUSTER
    numberOfClusters = 1
    abundanceFactors=[]
    numberClusterSpecificProteins=[]
    clusterSpecificProteinNames=[]
    scaleSameProteins=False
    
    #TRAJECTORIES
    """ We can simulate also trajectories instead of discrete clusters:
    for simplicity we CAN NOT specify those trajectories in the simulation.ini file,
    we can only set those trajectories ON, and then we have a hard-coded simple way of simulating those
    (all we can define are the number of proteins involved in those trajectories in the ini file)"""
    trajectory = None

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
            elif(str.startswith(line, "MinProtein=")):
                info = re.match(("MinProtein=(.*)"), line)
                self.MinProtein = int(info[1].rstrip("\n"))
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
                
            #CORRELATION VARIABLES 1.)
            elif(str.startswith(line, "randomProteinCorrelations=")):
                info = re.match(("randomProteinCorrelations=(.*)"), line.rstrip('\n'))
                info = info[1].rstrip("\n")
                self.randomProteinCorrelations = int(info)
                assert self.randomProteinCorrelations==1 or self.randomProteinCorrelations==0, "randomProteinCorrelations can only be 0,1 for true/false of random protein values (keep in mind, you ether set this OR give correlation specific distributions to sample from)"
            elif(str.startswith(line, "betaParameter=")):
                info = re.match(("betaParameter=(.*)"), line.rstrip('\n'))
                info = info[1].rstrip("\n")
                self.betaParameter = float(info)
                assert self.randomProteinCorrelations>0, "beta parameter must be > 0 (keep in mind, you ether set this OR give correlation specific distributions to sample from)"
            elif(str.startswith(line, "betaParameters=")):
                info = re.match(("betaParameters=(.*)"), line.rstrip('\n'))
                info = info[1].rstrip("\n")
                assert info[0] == "[", "First character must be '[' for following line in the ini file: " + line
                assert info[-1] == "]", "Last character must be ']' for following line in the ini file: " + line                
                self.addBetaParameters = ast.literal_eval(info)
                self.addBetaParameters = [float(element) for element in self.addBetaParameters]
            #CORRELATION VARIABLES 2.)
            elif(str.startswith(line, "proteinCorrelationDists=")): 
                #same for as ProteinLevels: [[mu, sd, #proteins]; ...] 
                # with as many Correlation Distributions as wanted
                info = re.match(("proteinCorrelationDists=\[(.*)\]"), line)
                info = str(info[1]).split(";")
                for element in info:
                    parts = re.match(("\s*\[\s*(-*\d+[\.\d+]*)\s*,\s*(\d+[\.\d+]*)\s*,\s*(\d+[\.\d+]*)\s*\]\s*"), element)
                    newCorrDist = self.ProteinCorrelation(float(parts[1]), float(parts[2]), int(float(parts[3])))
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
            elif(str.startswith(line, "trajectory=")): 
                trajectory = re.match(("trajectory=(.*)"), line)
                trajectory = str(trajectory[1])
                if( (trajectory == "linear") or (trajectory == "branche") or (trajectory == "cycle")):
                    self.trajectory = trajectory
                else:
                    print("Wrong trajectory stated:" + trajectory)
                    print("You must choose between: linear, branche, cycle!!\n")
                    exit()
            elif(str.startswith(line, "scaleSameProteins=")):
                scaleSameProteins = re.match(("scaleSameProteins=(.*)"), line)
                scaleSameProteins = int(scaleSameProteins[1])
                if(scaleSameProteins == 1):
                    self.scaleSameProteins = True

            line = file.readline()
        #we need exactly one abundance/ correlation/ percentage factor for every clusters
        self.numberOfClusters = len(self.cellPercentages)
        #set cellPercentages to 100 if no clusters exist
        if(self.cellPercentages == []):
            self.cellPercentages = [100]
            self.numberOfClusters = 1
            
        #CAUTION insufficient tests: we should test that ONLY 1.) or 2.) parameters are given, this is just a simple check
        # that when randomCorrelations is set we dont have all the other parameters set as well...
        if(self.randomProteinCorrelations==1):
            assert self.betaParameter > 0, "betaParameter must be above 0 if we want to make random correlations"
        if(self.randomProteinCorrelations==1):
            assert self.proteinCorrelationDists == [], "for random correlations proteinCorrelationDists does not have to be set"
            assert self.proteinCorrelationMean == 0.0, "for random correlations proteinCorrelationMean does not have to be set"
            assert self.proteinCorrelationStd == 0.4, "for random correlations proteinCorrelationStd does not have to be set"
            assert self.correlationSets == [], "for random correlations correlationSets does not have to be set"
            assert self.correlationFactors == [], "for random correlations correlationFactors does not have to be set"

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
        self.filePath = paramter_file

class SingleCellSimulation():

    #certain variables drawn from distributions (correlations) have to be stored to be accessed in benchmark
    """ BenchmarkIniFiles """
    benchmarkIniPath = ""     

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
        indexSlash = self.parameters.filePath.rfind("/")
        self.benchmarkIniPath = self.parameters.filePath[:indexSlash + 1] + "benchmarkIniFile" + self.parameters.filePath[indexSlash + 1:]

    #e.g., to write the proteins involved in clusters to benchmarkIniFile
    def __write_cluster_proteins(self, file_path, data):
        """
        Write a list to a file.
        :param file_path: The path to the file.
        :param data: The list of lists: protein list per cluster
        """
        try:
            # Open the file in write mode
            with open(file_path, 'w') as file:
                # Write each element of the list to the file
                for idx in range(self.parameters.numberOfClusters-1):
                    file.write(str(idx) + ':\t')
                    for protein in data[idx]:
                        file.write(str(protein) + '\t')
                    file.write('\n')
            print("Data written to", file_path)
        except Exception as e:
            print("Error writing to file:", str(e))
            
    #e.g., to write the proteinCorrelations to benchmarkIniFile (witha  cutoff stated in main)
    def __write_correlated_proteins(self, file_path, covariancematrix):
        """
        Write a list to a file.
        :param file_path: The path to the file.
        :param data: The list to write to the file.
        """
        try:
            # Open the file in write mode
            with open(file_path, 'w') as file:
                # Write each element of the list to the file
                print(file_path)
                for item in data:
                    print(item)

                    file.write(str(item) + '\n')
            print("Data written to", file_path)
        except Exception as e:
            print("Error writing to file:", str(e))
            
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
        #sample lbisize factors from log normal dist: 1.value = mean, 2.= std
        perturbFrame["factor"] = np.random.lognormal(self.parameters.libSize[0],self.parameters.libSize[1], len(perturbFrame))
        print(perturbFrame)
        perturbDict = dict(zip(perturbFrame["sample_id"], perturbFrame["factor"]))
        data["ab_count"] = data["ab_count"] / data["sample_id"].map(perturbDict)
        self._libsizeFactors = perturbDict

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
        std = 0.01
        if(setIdx is None):
            mean = self.parameters.proteinCorrelationMean * scalingFactor
            std = self.parameters.proteinCorrelationStd
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
    
    #method according to Lewandowski et. al 2009
    # better than second method (from post see below) since it also introduces neg. correlations
    # implementation based on matlab code from https://stats.stackexchange.com/questions/124538/how-to-generate-a-large-full-rank-random-correlation-matrix-with-some-strong-cor
    # 0.2 - 10 is a 'good' parameter range for strong to weak correlations for 60 proteins (manually observed)
    def __vineBeta(self, d, betaparam):
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
        minEigenVal = min(np.linalg.eigvals(covMat))
        if(minEigenVal <0):
            covMat = nearcorr(np.array(covMat), max_iterations=10000)

        return(covMat)
        
    #make a list of all the protein distributions: defined by mu/ size
    def __generateProteinDistributions(self):
            
        distributions = []
        protein = 0

        for proteinRange in self.parameters.ProteinLevels:
            
            #generate mu/ size for all proteins of the proteinRange distribution (sampled from distribution of means and the calculating size)
            dist = ProteinCountDistribution(proteinRange.number, proteinRange.mu, proteinRange.size, self.parameters.MinProtein)
            proteinCountVector = dist.distributionValues() #vector of mean values for all proteins for this proteinRange-distribution
             
            for i in range(proteinRange.number):       
                
                #get the very specific mu/ size for the protein i in this proteinRange-distribution
                mu_tmp = proteinCountVector[i]
                size_tmp = self.__calculate_size_from_mu(mu_tmp)
                         
                tmpDict = {'mu': mu_tmp, 'size': size_tmp}
                distributions.append(tmpDict)

                protein += 1

        return(distributions)
    
    #add cluster specific effects and convert to n, p
    def __generateClusterSpcecificProteinDistributions(self,baseDistributions, scaledProteinIdxs, clusterIdx):
            
        distributions = []
        protein = 0

        proteinsToScale = scaledProteinIdxs[clusterIdx]
        for dist in baseDistributions:
                            
                #get the very specific mu/ size for the protein i in this proteinRange-distribution
                mu_tmp = dist["mu"]
                size_tmp = dist["size"]
                         
                #scale mu by cluster specific scaling factor
                if(protein in proteinsToScale):
                    assert clusterIdx != 0, "The zero Cluster(baseline, not sclaed) can not have a scaling factor"
                    mu_tmp = mu_tmp * self.parameters.abundanceFactors[clusterIdx]

                n,p = convert_neg_binom_params(mu_tmp, size_tmp)
                tmpDict = {'loc': 0.0, 'n': n, 'p': p, 'type': NegBinomUnivariate}
                distributions.append(tmpDict)

                protein += 1

        return(distributions)
    
    def __sampleCellTimesProteinMatrix(self, percentage, dist, cov, clusterIdx):
        
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
        infNum = np.isinf(copulaResult).values.sum()
        if(infNum > 0):
            print("simulated data contains " + infNum + "INF values: check simulations! (e.g., wrong correlation parameters\n)")
            exit()
            
        copulaResult = copulaResult.melt(id_vars = ["sample_id"], var_name="ab_id", value_name="ab_count")
        copulaResult.insert(0, 'batch_id', 'batch')
        copulaResult.insert(0, 'ab_type', 'any_ab')

        #insert treatment/cluster column +
        #add cluster prefix to sample_id (since we have the same idx for all clusters)
        if(clusterIdx == 0):
            copulaResult.insert(0, 'cluster_id', 'control')
            copulaResult['sample_id'] = 'control_' + copulaResult['sample_id'].astype(str)
        else:
            copulaResult.insert(0, 'cluster_id', ('cluster_' + str(clusterIdx)))
            copulaResult['sample_id'] = 'cluster_' + str(clusterIdx) + "_" + copulaResult['sample_id'].astype(str)

        return(copulaResult)
    
    #1.) random correlaitons according to lewandowski        
    def __generateRandomCorrelations(self):
        covariancematrix = []
        
        betaParams = self.parameters.addBetaParameters
        betaParams.insert(0, self.parameters.betaParameter)
        assert(len(betaParams) == self.parameters.numberOfClusters)
        for idx in range(self.parameters.numberOfClusters):
            covariancematrix.append(self.__vineBeta(self.parameters.ProteinNumber, betaParams[idx]))

        return(covariancematrix)
    
    #correlations drawn from specific distributions
    def __generateSampledCorrelations(self):
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

        #if there r not corelation differences between cluster, generate ONLY ONE COV matrix
        #otherwise slight corr differences will result in perfect cluster seperation
        if not len(correlationProteinsSets):
            idx = 0
            covMatTmp = self.__generateCovarianceMatrix(correlationProteinsSets, scaledCorrelationIdxs, correlationScalingFactors, idx)
            for idx in range(self.parameters.numberOfClusters):
                covariancematrix.append(covMatTmp)
        else:
            for idx in range(self.parameters.numberOfClusters):
                covariancematrix.append(self.__generateCovarianceMatrix(
                    correlationProteinsSets, scaledCorrelationIdxs, correlationScalingFactors, idx))

        return(covariancematrix)
    
    def __generateClusterSpecificCovarianceMatrices(self):
        covariancematrix = []
        
        if(self.parameters.randomProteinCorrelations == 1):
            covariancematrix = self.__generateRandomCorrelations()
        else:
            covariancematrix = self.__generateSampledCorrelations()

        return(covariancematrix)
    
    def add_AB_to_name(self, list):
        # Iterate through each sublist
        for i in range(len(list)):
            # Iterate through each number in the sublist
            for j in range(len(list[i])):
                # Concatenate "AB" with the number and replace the original number
                list[i][j] = "AB" + str(list[i][j])
                
    def __sigmoidal_protein_activation(self, x):
        #sigmoidal function goes from -INF to INF
        #most of the change happens however between -4 and 4
        # we scale sigmoidal between 0-1 (this is then the factor for the protein: e.g., protein has a fold cahnge of 2, the actual count is then 
        # control + sigmoid(x)*fold-change * control )
        # scaling works with: (value - 0.5) * 8 to go from 0_1 range to -4_4
        
        scaled_x = (x - 0.5) * 8
        return(1 / (1 + np.exp(-scaled_x)))
    
    def __generateStartDict(self, overlap, pseudoTimeProteins):
        startDict = {}
        start = 0.0
        for p in pseudoTimeProteins:
            startDict[p] = start
            start += overlap
        return(startDict)
    
    def __pseudotime_protein_count(self, pseudoTimeProteins, baseDistributions, startDict, 
                                   proteinFactor, timeStepLength, cellFraction):
        #result data frame
        data = {
            'sample_id': [],
            'ab_id': [],
            'ab_count': []
        }
        pseudoTimeProteinCounts = pd.DataFrame(data)

        #follow pseudo-time and sample proteins
        cellnumber = 0
        for timeStep in range(1, int(self.parameters.CellNumber * cellFraction) + 1):
            timePoint = timeStep * timeStepLength
            for protein in pseudoTimeProteins:
                sigmoidalActivationFactor = 0
                if(timePoint>=startDict[protein]):
                    sigmoidalActivationFactor = self.__sigmoidal_protein_activation(timePoint-startDict[protein])

                #here sample the p[rotein count form a negative binomial (only take a single sample)
                proteinCountMean = baseDistributions[protein].get("mu") + sigmoidalActivationFactor * proteinFactor * baseDistributions[protein].get("mu")
                dist = ProteinCountDistribution( 1, mu = proteinCountMean, 
                                                size = baseDistributions[protein].get("size"), MinProtein = self.parameters.MinProtein)
                proteinCount = dist.distributionValues()[0]
                
                proteinCountDict = {'sample_id': 'sample_'+ str(cellnumber+1), 'ab_id': "AB"+str(protein), 'ab_count': proteinCount}
                proteinCountFrame = pd.DataFrame([proteinCountDict])
                pseudoTimeProteinCounts = pd.concat([pseudoTimeProteinCounts, proteinCountFrame], ignore_index=True)
            cellnumber += 1
        return(pseudoTimeProteinCounts)   
    
    def __pseudotime_protein_count_with_degradation(self, pseudoTimeProteins, baseDistributions, startDict, 
                                   proteinFactor, timeStepLength, cellFraction):
        #result data frame
        data = {
            'sample_id': [],
            'ab_id': [],
            'ab_count': []
        }
        pseudoTimeProteinCounts = pd.DataFrame(data)
        
        #total timeFrame length
        timeFrameLength = int(self.parameters.CellNumber * cellFraction) * timeStepLength

        #follow pseudo-time and sample proteins
        cellnumber = 0
        print(startDict)
        
        for timeStep in range(1, int(self.parameters.CellNumber * cellFraction) + 1):
            timePoint = timeStep * timeStepLength
            print("TIME: "+str(timePoint))
            print("_____")
            for protein in pseudoTimeProteins:
                print(str(protein) + ": ",end='')
                sigmoidalActivationFactor = 0
                #reset timePoint to fit into the infinite timeBox (imagine t = 0 is also t=max since we enter t=0 at the end our time frame)
                infiniteTimePoint = timePoint
                #check if we have to re-set the infinite timePoint
                
                #1.)if the startPoint of the protein is behind the timepoint (timeLength - 2), it must affect my current time point
                #and we need to temporaryly reset our infinite timeFrame
                #2 for activation + deactivation eachj of length 1
                
                #2.) if the current time point is actually between the proteinStart and timeFrameEnd then we r right now
                #within the activation/ deactivation range of the protein, even though this protein also spans
                #then infinite time point and we here have to treat it as a normal protein
                if( (timeFrameLength - 2) < (startDict[protein]) and 
                   not ( ((startDict[protein])<=timePoint) and (timePoint <=timeFrameLength ))  ):
                    infiniteTimePoint = timeFrameLength + timePoint
                    print("infinite", end = '')
                
                if(infiniteTimePoint >= startDict[protein] and infiniteTimePoint < (startDict[protein]+1.0) ):
                    sigmoidalActivationFactor = self.__sigmoidal_protein_activation(infiniteTimePoint-startDict[protein])
                    print(" up ",end='')

                elif(infiniteTimePoint >= (startDict[protein]+1.0) and infiniteTimePoint <= (startDict[protein]+2.0)):
                    #sigmoidal deactivation is basically 1 - activation
                    sigmoidalActivationFactor = (1 - self.__sigmoidal_protein_activation(infiniteTimePoint-(startDict[protein]+1)))
                    print(" down ",end='')

                proteinCountMean = baseDistributions[protein].get("mu") + sigmoidalActivationFactor * proteinFactor * baseDistributions[protein].get("mu")
                dist = ProteinCountDistribution( 1, mu = proteinCountMean, 
                                                size = baseDistributions[protein].get("size"), MinProtein = self.parameters.MinProtein)
                proteinCount = dist.distributionValues()[0]
                
                proteinCountDict = {'sample_id': 'sample_'+ str(cellnumber+1), 'ab_id': "AB"+str(protein), 'ab_count': proteinCount}
                proteinCountFrame = pd.DataFrame([proteinCountDict])
                pseudoTimeProteinCounts = pd.concat([pseudoTimeProteinCounts, proteinCountFrame], ignore_index=True)
            cellnumber += 1
        return(pseudoTimeProteinCounts)  
    
    """ sample all the proteins from a neg. binom. dist except the ones in
    pseudoTimePropteins, as those were sampled to follow the specific pseudoTime
    trajectory"""    
    def __sample_proteins_from_negBinom(self, baseDistributions, pseudoTimeProteins, cellFraction):
        #the result matrix
        proteinCountMatrix = None
        #for every proteins except the pseudotime ones
        #firstly store it in a df with one column per AB, then pivot it in the end
        for protIdx in range(self.parameters.ProteinNumber):
            if(not protIdx in pseudoTimeProteins):
                dist = ProteinCountDistribution( int(self.parameters.CellNumber * cellFraction), mu = baseDistributions[protIdx].get("mu"), 
                                                size = baseDistributions[protIdx].get("size"), MinProtein = self.parameters.MinProtein)
                proteinCountVector = dist.distributionValues()
                abName = "AB" + str(protIdx)
                #add it to a matrix as a new column
                if(proteinCountMatrix is None):
                    proteinCountMatrix = pd.DataFrame({abName : proteinCountVector}) 
                    proteinCountMatrix.index = ["sample_" + str(j+1) for j in range( int(self.parameters.CellNumber * cellFraction))]
                else:
                    proteinCountMatrix[abName] = proteinCountVector

        #convert dataframe into final format
        proteinCountMatrix = proteinCountMatrix.reset_index().rename(columns={ 'index' : 'sample_id'})
        proteinCountMatrix = proteinCountMatrix.melt(id_vars = ["sample_id"], var_name="ab_id", value_name="ab_count")
        return(proteinCountMatrix)             
                
    def __generateLinearGroundTruth(self, clusterSpecificProteinIdxs, overlap = 0.25,
                                    cellFraction = 1.0, baseDistributions = None):
        #generate pseudo-time function (one function of overlapping sigmoidals mapping time to locatio in protein space)
        
        #calculate time-steps 
        #pseudotime runs from 0 to sum of all sigmoidals per protein (overlapping)
        #overlap of 0.5 means next sigmoid start after previous one was activated to 50%
        numberPorteins = self.parameters.numberClusterSpecificProteins[0] #proteins in first cluster (we should only have ONE additional clsuter, numberClusterSpecificProteins does not contain control)
        #summed time are all the sigmoid parts before the overlap of the next plus what is missing in the end for the sigmoid
        summedPseudoTime = overlap * numberPorteins + (1.0 - overlap)
        #individual time steps from cell to cell
        timeStepLength = summedPseudoTime/ float(self.parameters.CellNumber * cellFraction)
        #both clusterSpecificProteinIdxs and self.parameters.abundanceFactors has control at 0 index
        pseudoTimeProteinList = clusterSpecificProteinIdxs[1] #list of proptein indices
        proteinFactor = self.parameters.abundanceFactors[0] #same factor for all proteins of same cluster: 0 is first cluster and not control
        
        print(pseudoTimeProteinList)
        #sample mean of all proteins
        if(baseDistributions is None):
            baseDistributions = self.__generateProteinDistributions() 
        # dicts map the INDEX to the MEAN/START-TIME
        startDict = self.__generateStartDict(overlap, pseudoTimeProteinList)
        
        #calculate actual protein count after every pseudotime-step
        pseudoTimeProteins = self.__pseudotime_protein_count(pseudoTimeProteinList, baseDistributions, startDict, 
                                                             proteinFactor, timeStepLength, cellFraction)
        pseudoTimeProteins.insert(0, 'ab_type', 'trajectoryProteins')
        
        #calculate all the other random proteins
        randomProteins = self.__sample_proteins_from_negBinom(baseDistributions, pseudoTimeProteinList, cellFraction)
        randomProteins.insert(0, 'ab_type', 'randomProteins')

        #combine data frames
        linearGroundTruth = pd.concat([pseudoTimeProteins, randomProteins], ignore_index = True)
        linearGroundTruth["ab_count"] = np.round(linearGroundTruth["ab_count"])
        linearGroundTruth.insert(0, 'cluster_id', 'no_cluster')
        linearGroundTruth.insert(0, 'batch_id', 'no_batch')
        print(linearGroundTruth)
        
        return(linearGroundTruth)
                    
    """ proteins accumulate until ONE cell state and are then subsequentyl degraded again..."""
    def __generateCyclicGroundTruthWithLibsize(self, clusterSpecificProteinIdxs):
        #basically its the same as doing two linear ground truths that connect
        #linear ground truth is actually sequntial sigmoidal activations, we can run it twice in
        #inverted order of proteins, and with less overlap of sigmoidals to create a wider circle
        
        #make sure we use the same base dist for both trajectories (forward & backwards)
        baseDistributions = self.__generateProteinDistributions() 
        forwardOrder = self.__generateLinearGroundTruth(clusterSpecificProteinIdxs, overlap = 0.75, 
                                                        cellFraction = 0.5, baseDistributions = baseDistributions)
        #reverse protein Order for cluster specific proteins (cluster marks to maximum point of circle connected to control)
        clusterSpecificProteinIdxs[1] = clusterSpecificProteinIdxs[1][::-1]
        reverseOrder = self.__generateLinearGroundTruth(clusterSpecificProteinIdxs, overlap = 0.75, 
                                                        cellFraction = 0.5, baseDistributions = baseDistributions)
        #sample_ids now exist twice, we need to change the number
        ids = reverseOrder['sample_id'].str.extract(r'_(\d+)').astype(int)
        print(ids)
        ids = ids + int(self.parameters.CellNumber * 0.5)
        reverseOrder["sample_id"] = "sample_" + ids.astype(str)
        print(reverseOrder)
        print(reverseOrder["sample_id"])

        circleData = pd.concat([forwardOrder, reverseOrder], ignore_index = True)
        return(circleData)
    
    
    """ proteins are one after the other acitvated and immeidately degraded.
        During activation of first protein the last one is also already degraded"""
    def __generateCyclicGroundTruth(self, clusterSpecificProteinIdxs, overlap = 0.5,
                                    cellFraction = 1.0, baseDistributions = None):
        #generate pseudo-time function (one function of overlapping sigmoidals mapping time to locatio in protein space)
        
        #calculate time-steps 
        #pseudotime runs from 0 to sum of all sigmoidals per protein (overlapping)
        #overlap of 0.5 means next sigmoid start after previous one was activated to 50%
        numberPorteins = self.parameters.numberClusterSpecificProteins[0] #proteins in first cluster (we should only have ONE additional clsuter, numberClusterSpecificProteins does not contain control)
        #summed time are all the sigmoid parts before the overlap of the next plus what is missing in the end for the sigmoid
        #the activation is of length 1 but the deactivation as well!!!
        summedPseudoTime = overlap * numberPorteins
        #individual time steps from cell to cell
        timeStepLength = summedPseudoTime/ int(self.parameters.CellNumber * cellFraction)
        #both clusterSpecificProteinIdxs and self.parameters.abundanceFactors has control at 0 index
        pseudoTimeProteinList = clusterSpecificProteinIdxs[1] #list of proptein indices
        proteinFactor = self.parameters.abundanceFactors[0] #same factor for all proteins of same cluster: 0 is first cluster and not control
        
        print(pseudoTimeProteinList)
        #sample mean of all proteins
        if(baseDistributions is None):
            baseDistributions = self.__generateProteinDistributions() 
        # dicts map the INDEX to the MEAN/START-TIME
        startDict = self.__generateStartDict(overlap, pseudoTimeProteinList)
        
        #calculate actual protein count after every pseudotime-step
        pseudoTimeProteins = self.__pseudotime_protein_count_with_degradation(pseudoTimeProteinList, baseDistributions, startDict, 
                                                             proteinFactor, timeStepLength, cellFraction)
        pseudoTimeProteins.insert(0, 'ab_type', 'trajectoryProteins')
        
        #calculate all the other random proteins
        randomProteins = self.__sample_proteins_from_negBinom(baseDistributions, pseudoTimeProteinList, cellFraction)
        randomProteins.insert(0, 'ab_type', 'randomProteins')

        #combine data frames
        linearGroundTruth = pd.concat([pseudoTimeProteins, randomProteins], ignore_index = True)
        linearGroundTruth["ab_count"] = np.round(linearGroundTruth["ab_count"])
        linearGroundTruth.insert(0, 'cluster_id', 'no_cluster')
        linearGroundTruth.insert(0, 'batch_id', 'no_batch')
        print(linearGroundTruth)
        
        return(linearGroundTruth)

#    def __generateBranchingGroundTruth(self, clusterSpecificProteinIdxs):

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
            if(self.parameters.scaleSameProteins == 0):
                #for seperate proteins sample the set per cluster
                for idx in range(self.parameters.numberOfClusters-1):
                    subset = random.sample(range(self.parameters.ProteinNumber), self.parameters.numberClusterSpecificProteins[idx])
                    clusterSpecificProteinIdxs.append(subset)
            else:
                #otherwise sample ONCE and every cluster has the same proteins with fold-change
                firstClusterIdx = 0
                subset = random.sample(range(self.parameters.ProteinNumber), self.parameters.numberClusterSpecificProteins[firstClusterIdx])
                for idx in range(self.parameters.numberOfClusters-1):
                    if(idx != self.parameters.numberClusterSpecificProteins[firstClusterIdx]):
                        print("If all cluster have THE SAME PROTEINS, all clusters MUST also ahve the SAME NUMBER OF FOLD-CHANGE PROTEINS \
                              in the variable <numberClusterSpecificProteins>")
                        exit()
                    clusterSpecificProteinIdxs.append(subset)
        
        self.parameters.clusterSpecificProteinNames = copy.deepcopy(clusterSpecificProteinIdxs)
        self.add_AB_to_name(self.parameters.clusterSpecificProteinNames) 
            
        #CLUSTER SPECIFIC PROTEIN-DISTRIBUTION-PARAMETRIZATION (scaled for proteins)
        #firstly generate distributions for every proteins (mu/ size) independant of clusters
        baseDistributions = self.__generateProteinDistributions() 
        #convert neg.binom to n,p and additionally add scaling factors for cluster
        #we need to insert 1.0 abudnace factor for baseline clusterIdx=0
        self.parameters.abundanceFactors.insert(0, 1.0) #base cluster has no sclaing of protein abudnance (so times factor 1.0)
        for idx in range(self.parameters.numberOfClusters):
            distributions.append(self.__generateClusterSpcecificProteinDistributions(baseDistributions, clusterSpecificProteinIdxs, idx))

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
        if(self.parameters.trajectory is None):
            self.groundTruthData = self.__generateGroundTruthWithCopula()
        else:
            clusterSpecificProteinIdxs = [[]]
            if(self.parameters.numberOfClusters > 1):
                #for seperate proteins sample the set per cluster
                for idx in range(self.parameters.numberOfClusters-1):
                    subset = random.sample(range(self.parameters.ProteinNumber), self.parameters.numberClusterSpecificProteins[idx])
                    clusterSpecificProteinIdxs.append(subset)
            if(self.parameters.trajectory == "linear"):
                assert(len(clusterSpecificProteinIdxs) == 2, "We need two clusters with 1 cluster with specific proteins for linear trajectory \
                    to model how proteins cahnge from control to cluster_1")
                self.groundTruthData = self.__generateLinearGroundTruth(clusterSpecificProteinIdxs) 
            elif(self.parameters.trajectory == "cycle"):
                assert(len(clusterSpecificProteinIdxs) == 2, "We need two clusters with 1 cluster with specific proteins for cycle trajectory \
                    to model how proteins cahnge from control to cluster_1 in cyclic way")
                self.groundTruthData = self.__generateCyclicGroundTruth(clusterSpecificProteinIdxs) 
            elif(self.parameters.trajectory == "branch"):
                assert(len(clusterSpecificProteinIdxs) == 4)
                self.groundTruthData = self.__generateBranchingGroundTruth(clusterSpecificProteinIdxs, "We need 4 clusters with 3 cluster with specific proteins for linear trajectory \
                    to model how proteins change linearly to cluster_1, and then rbanch into two seperate trajectories.")

        #not interesting at this point (initialy thought of AB duplicates for normalization)
        if(self.parameters.abDuplicates > 1):
            self.groundTruthData = self.__insert_ab_duplicates(self.groundTruthData)

        #scale counts for batch effects (NOT USED AT THIS POINT)
        perturbedData = self.groundTruthData.copy(deep=True)
        if(not (self.parameters.batchFactors is None)):
            perturbedData = self.__insert_batch_effect(perturbedData)
            self.__add_batch_effect_to_ground_truth(perturbedData)
        #add random noise to simulations
        tmp_simulatedData = perturbedData.copy(deep=True)
        if(self.parameters.noise > 0.0):
            tmp_simulatedData = self.__perturb(tmp_simulatedData)
            
        #scale counts for library size effects
        if(not (self.parameters.libSize[0]==1 and self.parameters.libSize[1]==1)):
            tmp_simulatedData = self.__insert_libsize_effect(tmp_simulatedData)

        #simulate AB binding efficiency
        #discard a fraction of proteinCounts as no AB has bound to them
        #tmp_simulatedData = self.__simulate_ab_binding(perturbedData)

        #simulate PCR amplification and sequencing
        #sampling with replacement to simulate PCR amplification as well as missing out on reads during washing/ sequencing
        #tmp_simulatedData = self.__simulate_sequencing_binding_2(perturbedData)

        #remove negative counts from data, and substitute with zero
        self.simulatedData = tmp_simulatedData
        self.simulatedData = self.simulatedData.applymap(replace_negatives)
        self.simulatedData["ab_count"] = self.simulatedData["ab_count"].round(decimals = 0)

        return(self.simulatedData)
    
    """ write all metadata to a file:
    e.g.: the proteins per cluster that have a log-fold change comapred to control population
    FORMAT: we have three columns: 
    VARIABLE (type of metadata, e.g., logFoldProtein)
    KEY: e.g.: the cluster for which this variable accounts, could also be another key later on which is protein specific, ...
    VALUE: the actual value of fold-change
    
    """
    def write_metadata(self):
        
        metadata = pd.DataFrame(columns=['VARIABLE', 'KEY', 'VALUE'])
        
        #ADD CLUSTER SPECIFIC PROTEIN FOLD CHANGES TO META DATA
        clusterNum = 0
        if(not len(self.parameters.clusterSpecificProteinNames)==0):
            for sublist in self.parameters.clusterSpecificProteinNames:
                # Iterate through each element in the sublist
                for element in sublist:
                    #create the new row
                    if(clusterNum == 0):
                        key = "control"
                    else:
                        key = "cluster_" + str(clusterNum)
                    new_row = {'VARIABLE': "LogFoldProtein", 'KEY': key, 'VALUE': str(element)}
                    metadata = metadata.append(new_row, ignore_index=True)
                clusterNum += 1
                
        #write all beta factors:
        # variable "LibsizeFactor"
        for factor, key in self._libsizeFactors.items():
            new_row = {'VARIABLE': "libsizeFactor", 'KEY': key, 'VALUE': str(factor)}
            metadata = metadata.append(new_row, ignore_index=True)
                
        #WRITE THE CLUSTER SPECIFIC PROTEINS
        metadata.to_csv(self.output_dir + "/" + self.parameters.simulationName + "_metadata.tsv", sep='\t', index = False)
        

    """ SAVE THE DATA """
    def save_data(self):
        #safe data
        printToTerminalOnce("\tSave Data\n")

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

        self.simulatedData.to_csv(self.output_dir + "/" + self.parameters.simulationName + "_SIMULATED.tsv", sep='\t', index = False)
        self.groundTruthData.to_csv(self.output_dir + "/" + self.parameters.simulationName + "_GROUNDTRUTH.tsv", sep='\t', index = False)
        
        printToTerminalOnce("\tData saved\n")
