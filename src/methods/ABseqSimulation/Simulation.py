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

#generates a ground truth matrix for protein counts
class Parameters():

    #incubation time is considered constant for now
    #antibody count is not considered for now: in sc IDseq constant AB concentration was used for all targets (0.1myg/ml over night)
    #biological variation matrix
    #technical variation matrix
    ProteinLevel = namedtuple('ProteinLevel', ['start', 'end', 'number'])

    """ SIMULATION Parameters """
    abBindingEfficiency = None
    seqAmplificationEfficiency = None
    #list of namedtuples ProteinLevel
    ProteinLevels = None
    size = None
    CellNumber = None
    abDuplicates = 1
    pcrCycles=1
    pcrCapture=1
    treatments=0
    treatmentVector = None
    batches=1
    batchFactors=None

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
            elif(str.startswith(line, "pcrCycles")):
                info = re.match(("pcrCycles=(.*)"), line)
                self.pcrCycles = float(info[1].rstrip("\n"))
            elif(str.startswith(line, "pcrCapture")):
                info = re.match(("pcrCapture=(.*)"), line)
                self.pcrCapture = float(info[1].rstrip("\n"))  
            elif(str.startswith(line, "treatments")):
                info = re.match(("treatments=(.*)"), line)
                self.treatments = float(info[1].rstrip("\n"))  
            elif(str.startswith(line, "treatmentVector")):
                info = re.match(("treatmentVector=(.*)"), line)
                info = str(info[1]).split(";")
                for elementVec in info:
                    treatmentAlteration = str(elementVec).split(",")
                    tmpTreatmentVec = None
                    for element in treatmentAlteration:
                        element = element.replace("[","") 
                        element = element.replace("]","") 
                        treatmentNum = float(element)
                        if(tmpTreatmentVec is not None):
                            tmpTreatmentVec.append(treatmentNum)
                        else:
                            tmpTreatmentVec = [treatmentNum]
                    if(self.treatmentVector is not None):
                        self.treatmentVector.append(tmpTreatmentVec)
                    else:
                        self.treatmentVector = [tmpTreatmentVec]
            elif(str.startswith(line, "batches")):
                            info = re.match(("batches=(.*)"), line)
                            self.batches = float(info[1].rstrip("\n"))  
            elif(str.startswith(line, "batchFactors")):
                            info = re.match(("batchFactors=(.*)"), line)
                            info = str(info[1]).split(",")
                            for batchNum in info:
                                num = float(batchNum)
                                if(self.batchFactors is not None):
                                    self.batchFactors.append(num)
                                else:
                                    self.batchFactors = [num]           
                                                         
    def __init__(self, paramter_file):
        self.__parseParameters(paramter_file)

class SingleCellSimulation():

    """ PARAMETERS """
    parameters = None
    ab_sampled = None
    output_dir = "bin/SIMULATIONS/"

    """ MODELLED DATA """
    groundTruthData = None
    #sample * AB matrix of simulated data
    simulatedData = None

    def __init__(self, parameters):
        self.parameters = parameters

    """ Generating GroundTruth of the Protein abundancies in all single cells"""
    def __generateGroundTruth(self, parameters):
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
        
        proteinCountMatrix = proteinCountMatrix.reset_index().rename(columns={ 'index' : 'sample_id'})
        proteinCountMatrix = proteinCountMatrix.melt(id_vars = ["sample_id"], var_name="ab_id", value_name="ab_count")
        self.groundTruthData = proteinCountMatrix

    """ Simulating the Detection of the GroudnTruth Protein Counts """
    def __simulate_ab_binding(self, data):
        number = int(self.parameters.abBindingEfficiency * len(data.index))
        tmp_simulatedData = data.sample(n=number, replace=False, random_state=1, weights = 'ab_count')
        return(tmp_simulatedData)

    def __generateUmiData(self, data):
        umiData = data.loc[data.index.repeat(data.ab_count)]
        umiData["umi_id"] = range(len(umiData.index))
        return(umiData)

    def __pcrAmplify(self, data):
    
        pcrNumber = int(self.parameters.pcrCapture * len(data.index))

        for i in range(int(self.parameters.pcrCycles)):
            print("PCR Cycle " + str(i))
            tmp_readsToPcrAmplify = data.sample(n=pcrNumber, replace = False, random_state=1)
            data = pd.concat([data, tmp_readsToPcrAmplify])
                
        #sample from all UMIs and remove umis that occur several times
        seqNumber = int(self.parameters.seqAmplificationEfficiency * len(data.index))
        data = data.sample(n=seqNumber, replace = False, random_state=1)

        data = data.drop_duplicates()

        return(data)
        
    def __simulate_sequencing_binding(self, data):
        #simulate UMIs
        #for every line simulate UMIs according to ab_count column
        print("Generate single reads per protein count with UMI")
        umiData = self.__generateUmiData(data)

        #PCR amplify those reads
        #sampling several times(for each PCR cycle) wo replacement and combine the old and sampled dataset
        # (previously simply sample w replacement and a sampling number > 1: however this neglects the idea that reads sampled in the first PCR
        # cycles are more likely to be overrepresented in the downrun)
        print("Generate PCR amplifications of reads")
        pcrData = self.__pcrAmplify(umiData)
        
        #concatenate again reads for same Protein in same cell after removing UMI column
        print("Concatenate all reads again to protein counts per cell")
        concatenatedData = pcrData.drop(columns=["umi_id"])
        concatenatedData["ab_count"] = concatenatedData.groupby(['sample_id', 'ab_id'])['sample_id'].transform('size')
        concatenatedData = concatenatedData.drop_duplicates()

        return(concatenatedData)

    """ MAIN FUNCTION: 
    1. generates ground truth & 
    2. simulates the protein count detection
    """
    def simulateData(self):

        #generate ground truth of the protein levels in each cell
        self.__generateGroundTruth(self.parameters)

        #simulate AB binding efficiency
        #discard a fraction of proteinCounts as no AB has bound to them
        tmp_simulatedData = self.__simulate_ab_binding(self.groundTruthData)

        #self.ab_sampled = tmp_simulatedData

        #simulate PCR amplification and sequencing
        #sampling with replacement to simulate PCR amplification as well as missing out on reads during washing/ sequencing
        tmp_simulatedData = self.__simulate_sequencing_binding(tmp_simulatedData)
        self.simulatedData = tmp_simulatedData

        return(self.simulatedData)

    """ SAVE THE DATA """
    def save_data(self, groundTruth = False):
        #safe data
        print("Safe Data")
        if(groundTruth):
            self.groundTruthData.to_csv(self.output_dir + "/GroundTruthData.tsv", sep='\t')
        else:
            self.simulatedData.to_csv(self.output_dir + "/SimulatedData.tsv", sep='\t')
            #self.ab_sampled.to_csv(self.output_dir + "/AbData.tsv", sep='\t')


