from dataclasses import replace
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

def LINE():
    return sys._getframe(1).f_lineno

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

    @dataclass
    class proteinCorrelation:
        prot1: int
        prot2: int
        positiveCorrelation: bool
        factor: int

    #incubation time is considered constant for now
    #antibody count is not considered for now: in sc IDseq constant AB concentration was used for all targets (0.1myg/ml over night)
    #biological variation matrix
    #technical variation matrix
    ProteinLevel = namedtuple('ProteinLevel', ['start', 'end', 'number'])

    """ SIMULATION Parameters """
    simulationName = ""
    abBindingEfficiency = None
    seqAmplificationEfficiency = None
    #list of namedtuples ProteinLevel
    ProteinLevels = None
    size = None
    CellNumber = None
    abDuplicates = 1
    abDuplicateRange = [1,1]
    abDuplicateDisturbance=0
    pcrCycles=1
    pcrCapture=1
    libSize=[1,1]
    treatmentVector = None
    diffExProteins = None
    batchFactors=None

    proteinCorrelations = []

    """ PARAMETERS
        rangeVector: a vector of quadruples(range, number, mean, std) for Abs of different distributions
        abDuplicates: to simulate the usage of several tags for one AB
    """

    def __parseParameters(self, paramFile):
        self.simulationName = os.path.basename(paramFile.rstrip('.ini'))
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
            elif(str.startswith(line, "diffExProteins")):
                info = re.match(("diffExProteins=(.*)"), line)
                info = str(info[1]).split(";")
                for elementVec in info:
                    diffExProteinIds = str(elementVec).split(",")
                    tmpProtIdVec = None
                    for element in diffExProteinIds:
                        element = element.replace("[","") 
                        element = element.replace("]","") 
                        protId = ("AB" + str(element))
                        if(tmpProtIdVec is not None):
                            tmpProtIdVec.append(protId)
                        else:
                            tmpProtIdVec = [protId]
                    if(self.diffExProteins is not None):
                        self.diffExProteins.append(tmpProtIdVec)
                    else:
                        self.diffExProteins = [tmpProtIdVec]
            elif(str.startswith(line, "batchFactors")):
                info = re.match(("batchFactors=(.*)"), line)
                info = str(info[1]).split(",")
                for batchNum in info:
                    num = float(batchNum)
                    if(self.batchFactors is not None):
                        self.batchFactors.append(num)
                    else:
                        self.batchFactors = [num]   
            elif(str.startswith(line, "libSize")):
                info = re.match(("libSize=\[(.*)\]"), line)
                info = str(info[1]).split(",")
                assert(len(info)==2)
                self.libSize[0]=float(info[0])
                self.libSize[1]=float(info[1])
            elif(str.startswith(line, "abDuplicateRange")):
                info = re.match(("abDuplicateRange=\[(.*)\]"), line)
                info = str(info[1]).split(",")
                assert(len(info)==2)
                self.abDuplicateRange[0]=float(info[0])
                self.abDuplicateRange[1]=float(info[1])
            elif(str.startswith(line, "proteinCorrelation=")):
                info = re.match(("proteinCorrelation=\[(.*)\]"), line)
                info = str(info[1]).split(",")
                for correlation in info:
                    p = re.compile('\[(\d*)([+-])(\d*)\]')
                    m = p.match(correlation)
                    posCorr = True
                    if(m[2] == "-"):
                        posCorr = False
                    cor = self.proteinCorrelation("AB"+(m[1]), "AB"+m[3], posCorr, 1)
                    self.proteinCorrelations.append(cor)
            elif(str.startswith(line, "proteinCorrelationFactors=")):
                info = re.match(("proteinCorrelationFactors=\[(.*)\]"), line)
                info = str(info[1]).split(",")
                i = 0
                for factor in info:
                    self.proteinCorrelations[i].factor = float(factor)
                    i+=1
            elif(str.startswith(line, "abDuplicateDisturbance=")):
                info = re.match(("abDuplicateDisturbance=(.*)"), line)
                self.abDuplicateDisturbance=float(info[1])

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
            dataTreatment["cluster_id"] = str(i)
            result.append(dataTreatment)
            i+=1
        concatedTreatments = pd.concat(result)
        self.groundTruthData = concatedTreatments

        return(concatedTreatments)

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
            dataBatch["batch_id"] = str(i)
            result.append(dataBatch)
            i+=1

        concatedBatches = pd.concat(result)
        return(concatedBatches)

    """ model libsize effect """
    def __insert_libsize_effect(self, data):
        print("Simulating different libsizes.")
        result = []
        sampleIdVector = (data["sample_id"].unique())
        for sampleId in sampleIdVector:
            n = random.uniform(self.parameters.libSize[0],self.parameters.libSize[1])
            dataSample = data[data["sample_id"] == sampleId]
            dataSample["ab_count"] = (n * dataSample["ab_count"])
            dataSample["ab_count"] = dataSample["ab_count"].round(decimals = 0)
            result.append(dataSample)

        concatedSamples = pd.concat(result)
        return(concatedSamples)

    """ model correlated proteins """
    def __insert_correlations_between_proteins(self, data):
        #sort data first according to sample
        #(when we change values for certain proteins we add values of another protein, the returned vector
        # might not be in the same order of sampels if not sorted first)
        data = data.sort_values(by=['sample_id'])

        for corr in self.parameters.proteinCorrelations:
            if(corr.positiveCorrelation):
                #origional negbinom distribution + factor times dependant protein count
                data.loc[data["ab_id"] == corr.prot2,"ab_count"] += np.array(data.loc[data["ab_id"] == corr.prot1,"ab_count"]) * corr.factor
            else:
                #origional negbinom dist - factor times dependant protein count (be aware: corr.factor is negative, so the whoel value is ADDED to old score)
                data.loc[data["ab_id"] == corr.prot2,"ab_count"] += (np.array(data.loc[data["ab_id"] == corr.prot1,"ab_count"]) * corr.factor)

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
    def __generateGroundTruth(self, parameters):
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

    def __add_batch_effect_to_ground_truth(self, newData):
        self.groundTruthData = pd.merge(self.groundTruthData,newData[['sample_id','ab_id','batch_id']],on=['sample_id','ab_id'], how='left')

    """ MAIN FUNCTION: 
    1. generates ground truth & 
    2. simulates the protein count detection
    """
    def simulateData(self):

        #generate ground truth of the protein levels in each cell
        self.__generateGroundTruth(self.parameters)
        #add additional features into ground truth data
        if(not (self.parameters.treatmentVector is None)):
            self.__insert_treatment_effect(self.groundTruthData)
        if(len(self.parameters.proteinCorrelations) > 0):
            self.__insert_correlations_between_proteins(self.groundTruthData)

        #add additional data pertubations
        perturbedData = self.groundTruthData.copy(deep=True)
        if(self.parameters.abDuplicates > 1):
            perturbedData = self.__insert_ab_duplicates(perturbedData)
        if(not (self.parameters.batchFactors is None)):
            perturbedData = self.__insert_batch_effect(perturbedData)
            self.__add_batch_effect_to_ground_truth(perturbedData)
        if(not (self.parameters.libSize[0]==1 and self.parameters.libSize[1]==1)):
            perturbedData = self.__insert_libsize_effect(perturbedData)

        

        #simulate AB binding efficiency
        #discard a fraction of proteinCounts as no AB has bound to them
        #deleted for now
        #tmp_simulatedData = self.__simulate_ab_binding(perturbedData)

        #self.ab_sampled = tmp_simulatedData

        #simulate PCR amplification and sequencing
        #sampling with replacement to simulate PCR amplification as well as missing out on reads during washing/ sequencing
        tmp_simulatedData = self.__simulate_sequencing_binding_2(perturbedData)

        self.simulatedData = tmp_simulatedData

        return(self.simulatedData)

    """ SAVE THE DATA """
    def save_data(self):
        #safe data
        print("Safe Data")

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

        
        self.simulatedData.to_csv(self.output_dir + "/" + self.parameters.simulationName + "_SIMULATED.tsv", sep='\t', index = False)
        self.groundTruthData.to_csv(self.output_dir + "/" + self.parameters.simulationName + "_GROUNDTRUTH.tsv", sep='\t', index = False)
