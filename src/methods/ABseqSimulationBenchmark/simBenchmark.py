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
import subprocess
import os
import shutil
import sys
sys.path.append('./src/methods/ToolBox')
from functions import *

class Parameters():

    iniFile = ""
    normMethods=[]
    batchEffect=False
    datasets = ""

    def __parseParameters(self, paramFile):
        self.iniFile = paramFile
        file = open(paramFile, "r")
        line = file.readline()
        while line:
            line = file.readline()
            if(str.startswith(line, "normMethods")):
                info = re.match(("normMethods=(.*)"), line)
                info = str(info[1]).split(",")
                for element in info:
                    self.normMethods.append(element)
            elif(str.startswith(line, "batchEffect")):
                info = re.match(("batchEffect=(.*)"), line)
                if(int(info[1]) == 1):
                    self.batchEffect = True

    def __parseDatasetDir(self):
        settingsFile = "./settings.ini"
        file = open(settingsFile, "r")
        line = file.readline()
        while line:
            line = file.readline()
            if(str.startswith(line, "datasets")):
                info = re.match(("datasets=(.*)"), line)
                self.datasets = str(info[1])

        if(self.datasets == ""):
            print("WARNING: NO DIRECTORY FOR NORMALIZATION FILES GIVEN")

    def __init__(self, paramFile):
        self.__parseParameters(paramFile)
        self.__parseDatasetDir()

class Benchmark():

    parameters=None

    def __init__(self, parameters):
        self.parameters = parameters

    def __generateNormalizationIni(self, normOriginFilePath):
        iniFile = removesuffix(normOriginFilePath,'.tsv') + ".ini"
        
        if os.path.exists(iniFile):
            os.remove(iniFile)
        iniStream = open(iniFile, "a")

        iniStream.write("sample_id=sample_id\n")
        iniStream.write("ab_id=ab_id\n")
        iniStream.write("ab_count=ab_count\n")
        iniStream.write("batch_id=batch_id\n")
        iniStream.write("cluster_id=cluster_id\n")
        iniStream.write("ab_type=ab_type\n")

        iniStream.write("removeBatchEffect=")
        if(self.parameters.batchEffect == True):
            iniStream.write(str(1))
        else:
            iniStream.write(str(0))

        iniStream.write("\n")

        iniStream.close()

    """ RUN a complete benchmark from data simulation up to Normalization evaluation """
    def run(self):
        #activate environment
        subprocess.run(["./source", "scRNAseq/bin/activate"], shell=True)

        #call simulations
        subprocess.run(["python3 ./src/methods/ABseqSimulation/main.py " + self.parameters.iniFile], shell = True, check = True)
        
        #copy simulations into normnalization folder:
        #from bin/SIMMULATIONS to datasets/
        simulationFilePath = "./bin/SIMULATIONS/"

        simulationName = os.path.basename(removesuffix(self.parameters.iniFile, '.ini'))
        simulatedFileName = simulationName + "_SIMULATED.tsv"
        
        normOriginFilePath = self.parameters.datasets + "/" + simulationName + ".tsv"
        simulationResultFilePath = simulationFilePath + simulatedFileName

        commandString = "cp " + simulationResultFilePath + " " + normOriginFilePath
        subprocess.run([commandString], shell = True, check = True)
        
        #generate a ini file next to it
        self.__generateNormalizationIni(normOriginFilePath)

        #run normalizations
        if( not self.parameters.normMethods ):
            print("No normalization methods given!")
            return
        for norm in self.parameters.normMethods:
            if(norm == "GRAPH"):
                commandString = "python3 src/methods/GraphNormalization/main.py " + normOriginFilePath
                try:
                    subprocess.run([commandString], shell = True, check = True)
                except:
                    print("ERROR: Normalization method " + norm + " failed")
            else:
                commandString = "Rscript --quiet ./src/normalization/NormalizationScript.R " + norm + " " + simulationName + ".tsv"
                try:
                    subprocess.run([commandString], shell = True, check = True)
                except:
                    print("ERROR: Normalization method " + norm + " failed")
        #move also ground truth into normalization folder
        groundTruthName = os.path.basename(removesuffix(self.parameters.iniFile, '.ini')) + "_GROUNDTRUTH.tsv"
        groundTruthResultFilePath = simulationFilePath + groundTruthName
        commandString = "cp " + groundTruthResultFilePath + " ./bin/NORMALIZED_DATASETS/" + simulationName + "/" + groundTruthName
        subprocess.run([commandString], shell = True, check = True)
        #move also simulated file into normalization folder
        simulatedName = os.path.basename(removesuffix(self.parameters.iniFile, '.ini')) + "_SIMULATED.tsv"
        simulatedResultFilePath = simulationFilePath + simulatedName
        commandString = "cp " + simulatedResultFilePath + " ./bin/NORMALIZED_DATASETS/" + simulationName + "/" + simulatedName
        subprocess.run([commandString], shell = True, check = True)

        try:
            #remove the ab_count with ab_count_normalized to still run benchmark on it
            #from ground truth
            groundTruthFile = ("./bin/NORMALIZED_DATASETS/" + simulationName + "/" + groundTruthName)
            groundTruthDataStream = open(groundTruthFile, "rt")
            dataGroundTruth = groundTruthDataStream.read()
            dataGroundTruth = dataGroundTruth.replace('ab_count', 'ab_count_normalized')
            groundTruthDataStream.close()
            groundTruthDataStream = open(groundTruthFile, "wt")
            groundTruthDataStream.write(dataGroundTruth)
            groundTruthDataStream.close()
            #from simulated file
            simulatedFile = ("./bin/NORMALIZED_DATASETS/" + simulationName + "/" + simulatedName)
            simulatedDataStream = open(simulatedFile, "rt")
            dataSimulated = simulatedDataStream.read()
            dataSimulated = dataSimulated.replace('ab_count', 'ab_count_normalized')
            simulatedDataStream.close()
            simulatedDataStream = open(simulatedFile, "wt")
            simulatedDataStream.write(dataSimulated)
            simulatedDataStream.close()
        except:
            print("ERROR for moving simulated and groundTruth data")

        #run normalization benchmark
        normResultFilePath = "./bin/NORMALIZED_DATASETS/" + simulationName
        benchmarkCommand = "python3 src/benchmark/normBenchmark.py --groundtruth --iniFile " + self.parameters.iniFile + " " + normResultFilePath
        subprocess.run([benchmarkCommand], shell = True, check = True)

    #every run of a simulation -> normalizaitons -> benchmark results in a folder in bin/BENCHMARK
    #all those folers are beeing stored in a foler called SIMULATIONS_dateTime to keep better track
    def moveIntoOneFolder(self, newSimulationDir):
        simulationName = os.path.basename(removesuffix(self.parameters.iniFile, '.ini'))
        if(os.path.exists("./bin/BENCHMARKED_DATASETS/" + simulationName)):
            shutil.move("./bin/BENCHMARKED_DATASETS/" + simulationName, newSimulationDir)

