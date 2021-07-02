from dataclasses import replace
from scipy.stats import nbinom
import re
import numpy as np
from collections import namedtuple
import random
import pandas as pd

#data from paper
    # 25 cycles of PCR

    #UMI DISTRIBUTION
    #size 3.83 && mu 4.68

    #PROTEIN DISTRIBUTION
    #most Abcounts are between 100 and 500, with a few beeing greater, and IGTB1 beeing around 23000 // 
    # size parameter has a median of 3.5 beeing mostly around 2.5 and 4

#introduce a biological variance matrix: e.g. neg. binomial distributed
#introduce a technical variance matrix: AB amount differences, annealing differences

class ProteinCountDistribution():
    abCountVector = None

    def __convert_params(self, mu, size):
        p = size/(size+mu)
        n = mu*p/(1.0 - p)
        return n, p

    def __init__(self, number, mu, size):
        #calculate n and p parameters
        n, p = self.__convert_params(mu, size)
        self.abCountVector = nbinom.rvs(n ,p, size = number)

    def distributionValues(self):
        return(self.abCountVector)


class Parameters():

    #incubation time is considered constant for now
    #antibody count is not considered for now: in sc IDseq constant AB concentration was used for all targets (0.1myg/ml over night)
    #biological variation matrix
    #technical variation matrix
    ProteinLevel = namedtuple('ProteinLevel', ['start', 'end', 'number'])

    """ MODELLED PARAMETERS """
    groundTruthData = None
    """ SIMULATION Parameters """
    abBindingEfficiency = None
    seqAmplificationEfficiency = None
    #list of namedtuples ProteinLevel
    ProteinLevels = None
    size = None
    CellNumber = None
    abDuplicates = 1

    """ PARAMETERS
        rangeVector: a vector of quadruples(range, number, mean, std) for Abs of different distributions
        abDuplicates: to simulate the usage of several tags for one AB
    """

    def __parseParameters(self, paramFile):
        file = open(paramFile, "r")
        line = file.readline()
        while line:
            line = file.readline()
            if(str.startswith(line, "ProteinLevels")):
                info = re.match(("ProteinLevels=\[(.*)\]"), line)
                info = str(info[1]).split(";")
                for element in info:
                    parts = re.match(("\s*\[\s*\[\s*([0-9]*)\s*,\s*([0-9]*)\s*\]\s*,\s*([0-9]*)\s*\]\s*"), element)
                    newProteinLevels = self.ProteinLevel(int(parts[1]), int(parts[2]), int(parts[3]))
                    if(self.ProteinLevels is not None):
                        self.ProteinLevels.append(newProteinLevels)
                    else:
                        self.ProteinLevels = [newProteinLevels]
            elif(str.startswith(line, "size")):
                info = re.match(("size=(.*)"), line)
                self.size = float(info[1].rstrip("\n"))
            elif(str.startswith(line, "CellNumber")):
                info = re.match(("CellNumber=(.*)"), line)
                self.CellNumber = int(info[1].rstrip("\n"))
            elif(str.startswith(line, "abDuplicates")):
                info = re.match(("abDuplicates=(.*)"), line)
                self.abDuplicates = int(info[1].rstrip("\n"))
            elif(str.startswith(line, "abBindingEfficiency")):
                info = re.match(("abBindingEfficiency=(.*)"), line)
                self.abBindingEfficiency = float(info[1].rstrip("\n"))
            elif(str.startswith(line, "seqAmplificationEfficiency")):
                info = re.match(("seqAmplificationEfficiency=(.*)"), line)
                self.seqAmplificationEfficiency = float(info[1].rstrip("\n"))

    def __simulateProteinCountMatrix(self, abDuplicates = 1):
        proteinCount = 1
        proteinCountMatrix = None
        for proteinRange in self.ProteinLevels:
            for i in range(proteinRange.number):
                #for every protein simulate a neg.binom distibuted number for every cell
                mu = random.randrange(proteinRange.start, proteinRange.end)
                #TO DO: size variance of neg.binom distribution is simpy a random value +- 1
                size_variance = random.randrange(-1, 1)
                size = self.size + size_variance
                dist = ProteinCountDistribution(self.CellNumber, mu, size)
                proteinCountVector = dist.distributionValues()
                abName = "AB" + str(proteinCount)
                #add it to a matrix as a new column
                if(proteinCountMatrix is None):
                    proteinCountMatrix = pd.DataFrame({abName : proteinCountVector}) 
                    proteinCountMatrix.index = ["sample_" + str(j+1) for j in range(self.CellNumber)]
                else:
                    proteinCountMatrix[abName] = proteinCountVector

                proteinCount +=1
        
        proteinCountMatrix = proteinCountMatrix.reset_index().rename(columns={ 'index' : 'sample_id'})
        proteinCountMatrix = proteinCountMatrix.melt(id_vars = ["sample_id"], var_name="ab_id", value_name="ab_count")
        self.groundTruthData = proteinCountMatrix

    def __init__(self, paramter_file):
        self.__parseParameters(paramter_file)
        self.__simulateProteinCountMatrix()

class SingleCellSimulation():

    parameters = None
    #sample * AB matrix of simulated data
    simulatedData = None
    ab_sampled = None
    output_dir = "bin/SIMULATIONS/"

    def __init__(self, parameters):
        self.parameters = parameters

    def __simulate_ab_binding(self, data):
        number = int(self.parameters.abBindingEfficiency * len(data.index))
        tmp_simulatedData = data.sample(n=number, replace=False, random_state=1, weights = 'ab_count')
        return(tmp_simulatedData)

    def __generateUmiData(self, data):
        umiData = data.loc[data.index.repeat(data.ab_count)]
        umiData["umi_id"] = range(len(umiData.index))
        return(umiData)

    def __simulate_sequencing_binding(self, data):
        #simulate UMIs
        #for every line simulate UMIs according to ab_count column
        umiData = self.__generateUmiData(data)

        number = int(self.parameters.seqAmplificationEfficiency * len(umiData.index))
        tmp_simulatedData = umiData.sample(n=number, replace = True, random_state=1)
        
        tmp_simulatedData["ab_count"] = tmp_simulatedData.groupby(['sample_id','ab_id'])['umi_id'].transform('size')
        tmp_simulatedData = tmp_simulatedData.drop(columns=["umi_id"])
        tmp_simulatedData = tmp_simulatedData.drop_duplicates()

        return(tmp_simulatedData)

    def simulateData(self):
        #simulate AB binding efficiency
        #discard a fraction of proteinCounts as no AB has bound them
        tmp_simulatedData = self.__simulate_ab_binding(self.parameters.groundTruthData)

        self.ab_sampled = tmp_simulatedData

        #simulate PCR amplification and sequencing
        #sampling with replacement to simulate PCR amplification as well as missing out on reads during washing/ sequencing
        tmp_simulatedData = self.__simulate_sequencing_binding(tmp_simulatedData)
        self.simulatedData = tmp_simulatedData

        return(self.simulatedData)

    def save_data(self, groundTruth = False):
        #safe data
        if(groundTruth):
            self.parameters.groundTruthData.to_csv(self.output_dir + "/GroundTruthData.tsv", sep='\t')
        else:
            self.simulatedData.to_csv(self.output_dir + "/SimulatedData.tsv", sep='\t')
            self.ab_sampled.to_csv(self.output_dir + "/AbData.tsv", sep='\t')


