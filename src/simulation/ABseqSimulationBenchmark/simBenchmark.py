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
            if(str.startswith(line, "normMethods")):
                info = re.match(("normMethods=(.*)"), line)
                info = str(info[1]).split(",")
                for element in info:
                    self.normMethods.append(element)
            elif(str.startswith(line, "batchEffect")):
                info = re.match(("batchEffect=(.*)"), line)
                if(int(info[1]) == 1):
                    self.batchEffect = True
            line = file.readline()

    def __parseDatasetDir(self):
        settingsFile = "./settings.ini"
        file = open(settingsFile, "r")
        line = file.readline()
        while line:
            if(str.startswith(line, "datasets")):
                info = re.match(("datasets=(.*)"), line)
                self.datasets = str(info[1])
            line = file.readline()

        if(self.datasets == ""):
            printToTerminalOnce("WARNING: NO DIRECTORY FOR NORMALIZATION FILES GIVEN")

    def __init__(self, paramFile):
        self.__parseParameters(paramFile)
        self.__parseDatasetDir()

class Benchmark():

    parameters=None

    def __init__(self, parameters, stdoutFile = "", noExplicitelySetThreads = False, keepData = False, benchmarkOnly = False, knnOverlap = 20, knnMetric = "manhattan"):
        self.parameters = parameters
        self.stdout = stdoutFile
        self.noExplicitelySetThreads = noExplicitelySetThreads
        self.keepData = keepData
        self.benchmarkOnly = benchmarkOnly
        self.knnOverlap = knnOverlap
        self.knnMetric = knnMetric

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
        subprocess.run(["source scRNAseq/bin/activate"], shell=True)

        #FOLDER/ FILE DEFINITION FOR SIMULATION
        simulationFilePath = "./bin/SIMULATIONS/"
        # ini files that simply maps the cooumn names from data.table to the unsed parameters in normalization
        simulationName = os.path.basename(removesuffix(self.parameters.iniFile, '.ini'))
        simulatedFileName = simulationName + "_SIMULATED.tsv"

        normOriginFilePath = self.parameters.datasets + "/" + simulationName + ".tsv"
        simulationResultFilePath = simulationFilePath + simulatedFileName

        folder_path = ("./bin/NORMALIZED_DATASETS/" + simulationName)

        #benchmark files
        groundTruthName = os.path.basename(removesuffix(self.parameters.iniFile, '.ini')) + "_GROUNDTRUTH.tsv"
        simulatedName = os.path.basename(removesuffix(self.parameters.iniFile, '.ini')) + "_SIMULATED.tsv"
        newGroundTruthName = "Groundtruth.tsv"
        newSimulatedName = "Simulation.tsv"
        groundTruthFile = ("./bin/NORMALIZED_DATASETS/" + simulationName + "/" + newGroundTruthName)
        simulatedFile = ("./bin/NORMALIZED_DATASETS/" + simulationName + "/" + newSimulatedName)
        metaDataFile = ("./bin/NORMALIZED_DATASETS/" + simulationName + "/metaData.tsv")

        if(self.benchmarkOnly):
            printToTerminalOnce("SKIP SIMULATIONS; go stright to benchmark: make sure INI files (check path in pipeline ini) and normalized folders per ini exist (in bin/NORMALIZED_DATASETS)")
        else:
            #1.) call simulations
            #__________________
            printToTerminalOnce("RUN SIMULATION")

            if(self.noExplicitelySetThreads):
                subprocess.run(["python3 ./src/simulation/ABseqSimulation/main.py " + self.parameters.iniFile + " --stdout " + self.stdout], shell = True, check = True)
            else:
                subprocess.run(["OMP_NUM_THREADS=1 USE_SIMPLE_THREADED_LEVEL3=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 ./src/simulation/ABseqSimulation/main.py " + self.parameters.iniFile + " --t 1 --stdout " + self.stdout], shell = True, check = True)

            #copy simulations into normnalization folder:
            ##################################
            #from bin/SIMMULATIONS to datasets/

            #copy simulated file in ./bin/SIMULATIONS/file_SIMULATWED.tsv into the folder for files to normalize
            #this file is stated in github dir ini file: for us 2_preprocessed
            commandString = "cp " + simulationResultFilePath + " " + normOriginFilePath
            subprocess.run([commandString], shell = True, check = True)

            #generate a ini file next to it
            self.__generateNormalizationIni(normOriginFilePath)

            #2.) run normalizations
            #__________________
            printToTerminalOnce("RUN NORMALIZATION")
            if( not self.parameters.normMethods ):
                printToTerminalOnce("No normalization methods given!")
                return

            #delete normalization folder if exists (since a previous mthod might have used other norm methods, with other cell numbers etc which might
            # lead to downstream errors)
            if os.path.exists(folder_path):
                shutil.rmtree(folder_path, ignore_errors=True)

            for norm in self.parameters.normMethods:
                if(norm == "GRAPH"):
                    commandString = ""
                    if(self.noExplicitelySetThreads):
                        commandString = "python3 src/methods/GraphNormalization/main.py " + normOriginFilePath + " " + folder_path + " " + str(0.5) + " >> " + self.stdout + " 2>&1"
                    else:
                        commandString = "OMP_NUM_THREADS=1 USE_SIMPLE_THREADED_LEVEL3=1 python3 src/methods/GraphNormalization/main.py " + normOriginFilePath + " " + folder_path + " " + str(0.5) + " >> " + self.stdout + " 2>&1"
                    try:
                        printToTerminalOnce("RUNNING: " + commandString)
                        subprocess.run([commandString], shell = True, check = True)
                    except Exception as e: 
                        print(e)
                        printToTerminalOnce("ERROR: Normalization method " + norm + " failed for " + self.parameters.iniFile)
                elif(norm == "SCVI"):
                    #for simplicity run this without threads (did not even try but had issues on LINUX for other methods before)
                    commandString = "OMP_NUM_THREADS=1 USE_SIMPLE_THREADED_LEVEL3=1 python3 src/methods/scVINormalization/main.py " + normOriginFilePath + " " + folder_path + " >> " + self.stdout + " 2>&1"
                    try:
                        printToTerminalOnce("RUNNING: " + commandString)
                        subprocess.run([commandString], shell = True, check = True)
                    except Exception as e: 
                        print(e)
                        printToTerminalOnce("ERROR: Normalization method " + norm + " failed for " + self.parameters.iniFile)
                elif(norm == "SANITY"):
                    #for simplicity run this without threads (did not even try but had issues on LINUX for other methods before)
                    commandString = "OMP_NUM_THREADS=1 USE_SIMPLE_THREADED_LEVEL3=1 python3 src/methods/SanityNormalization/main.py " + normOriginFilePath + " " + folder_path + " >> " + self.stdout + " 2>&1"
                    try:
                        printToTerminalOnce("RUNNING: " + commandString)
                        subprocess.run([commandString], shell = True, check = True)
                    except Exception as e: 
                        print(e)
                        printToTerminalOnce("ERROR: Normalization method " + norm + " failed for " + self.parameters.iniFile)            
                else:
                    commandString = ""
                    if(self.noExplicitelySetThreads):
                        commandString = "Rscript --quiet ./src/normalization/NormalizationScript.R " + norm + " " + simulationName + ".tsv >> " + self.stdout + " 2>&1"
                    else:
                        commandString = "OMP_NUM_THREADS=1 USE_SIMPLE_THREADED_LEVEL3=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --quiet ./src/normalization/NormalizationScript.R " + norm + " " + simulationName + ".tsv >> " + self.stdout + " 2>&1"
                    try:
                        printToTerminalOnce("RUNNING: " + commandString)
                        subprocess.run([commandString], shell = True, check = True)
                    except Exception as e: 
                        print(e)
                        printToTerminalOnce("ERROR: Normalization method " + norm + " failed for " + self.parameters.iniFile)

            try:

                #move also ground truth into normalization folder
                groundTruthResultFilePath = simulationFilePath + groundTruthName
                commandString = "cp " + groundTruthResultFilePath + " " + groundTruthFile
                subprocess.run([commandString], shell = True, check = True)
                #move also simulated file into normalization folder
                simulatedResultFilePath = simulationFilePath + simulatedName
                commandString = "cp " + simulatedResultFilePath + " " + simulatedFile
                subprocess.run([commandString], shell = True, check = True)
                
                #move meta-data into normalization folder
                metaDataFilePath = simulationFilePath + "/" + simulationName + "_metadata.tsv"
                commandString = "cp " + metaDataFilePath + " " + metaDataFile
                subprocess.run([commandString], shell = True, check = True)                

                #remove the ab_count with ab_count_normalized to still run benchmark on it
                #from ground truth
                groundTruthDataStream = open(groundTruthFile, "rt")
                dataGroundTruth = groundTruthDataStream.read()
                dataGroundTruth = dataGroundTruth.replace('ab_count', 'ab_count_normalized')
                groundTruthDataStream.close()
                groundTruthDataStream = open(groundTruthFile, "wt")
                groundTruthDataStream.write(dataGroundTruth)
                groundTruthDataStream.close()
                #from simulated file
                simulatedDataStream = open(simulatedFile, "rt")
                dataSimulated = simulatedDataStream.read()
                dataSimulated = dataSimulated.replace('ab_count', 'ab_count_normalized')
                simulatedDataStream.close()
                simulatedDataStream = open(simulatedFile, "wt")
                simulatedDataStream.write(dataSimulated)
                simulatedDataStream.close()
            except Exception as e: 
                print(e)
                printToTerminalOnce("ERROR for moving simulated and groundTruth data for " + self.parameters.iniFile + "\n")

        #run normalization benchmark
        try:
            normResultFilePath = "./bin/NORMALIZED_DATASETS/" + simulationName
            benchmarkCommand = ""
            if(self.noExplicitelySetThreads):
                benchmarkCommand = "python3 src/benchmark/main.py --groundtruth --iniFile " + self.parameters.iniFile + " --stdout " + self.stdout + " --t -1" + " --knn " + str(self.knnOverlap)  + " --metric " + self.knnMetric + " " + normResultFilePath
            else:
                benchmarkCommand = "OMP_NUM_THREADS=1 USE_SIMPLE_THREADED_LEVEL3=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 src/benchmark/main.py --groundtruth --iniFile " + self.parameters.iniFile + " --stdout " + self.stdout + " --t 1" + " --knn " + str(self.knnOverlap)  + " --metric " + self.knnMetric + " " + normResultFilePath
            printToTerminalOnce("Running benchmark with: " + benchmarkCommand + "\n")
            subprocess.run([benchmarkCommand], shell = True, check = True)
            
            #copy the raw data into the final folder
            finalFolderName = "./bin/BENCHMARKED_DATASETS/" + simulationName
            commandString = "cp -r " + normResultFilePath + " " + finalFolderName
            subprocess.run([commandString], shell = True, check = True)

            if(not self.keepData):
                #delete folder of normalized data
                subprocess.run("rm -r " + folder_path, shell = True, check = True)
                #delete simulated/ groundtruth data
                subprocess.run("rm -r " + simulationResultFilePath, shell = True, check = True)
                subprocess.run("rm -r " + groundTruthResultFilePath, shell = True, check = True)
                #delete the ini/tsv file copied into normalization folder
                iniFileToDelte = removesuffix(normOriginFilePath,'.tsv') + ".ini"
                subprocess.run("rm -r " + normOriginFilePath, shell = True, check = True)
                subprocess.run("rm -r " + iniFileToDelte, shell = True, check = True)
        except Exception as e: 
            print(e)
            printToTerminalOnce("ERROR Benchmarking of normalized files failed for: " + self.parameters.iniFile)

    #every run of a simulation -> normalizaitons -> benchmark results in a folder in bin/BENCHMARK
    #all those folers are beeing stored in a foler called SIMULATIONS_dateTime to keep better track
    def moveIntoOneFolder(self, newSimulationDir, dublicateNum):
        simulationName = str(os.path.basename(removesuffix(self.parameters.iniFile, '.ini')))
        if(os.path.exists("./bin/BENCHMARKED_DATASETS/" + simulationName)):
            if(os.path.exists(newSimulationDir + "/" + simulationName)):
                shutil.rmtree(newSimulationDir + "/" + simulationName)
            shutil.move("./bin/BENCHMARKED_DATASETS/" + simulationName, newSimulationDir)
            os.rename(newSimulationDir + "/" + simulationName, newSimulationDir + "/" + simulationName + "_" + str(dublicateNum))

        #copy the ini file 
        shutil.copy(self.parameters.iniFile, newSimulationDir + "/" + simulationName + "_" + str(dublicateNum))

    def combine_files(self, newSimulationDir, resultDir, fileName, dublicationNum):
        newFile = open(resultDir + fileName, 'a')  
        filePath = os.path.basename(removesuffix(self.parameters.iniFile, '.ini')) + "_" + str(dublicationNum)

        simulationName = newSimulationDir + "/" + filePath + "/Results/"
        #for loop to go through all files and get the ONE file with the results
        for oldFile in os.listdir(simulationName):
            if oldFile.endswith(fileName):                
                fileStream = open(simulationName + oldFile,'r')
                line = fileStream.readline()
                lineNum = 0
                while line:
                    if(lineNum==0 and os.stat(resultDir + fileName).st_size==0):
                        newFile.write("Simulation_Identifier" + "\t" + line)
                    elif(lineNum!=0):
                        fileIdentifier = re.findall(r'\d+_\d+', filePath)[0]
                        newFile.write(fileIdentifier + "\t" + line)
                    lineNum = lineNum + 1
                    line = fileStream.readline()
                fileStream.close()
                break
        newFile.close()

    #looks up the result files for all tests (Spearman/ Classification/ etv.) and writes the results
    #of the different Simulations into one common file
    def copyResultsIntoOneFile(self, newSimulationDir, dublicationNum):
        resultDir = newSimulationDir + "/Results/"
        if not os.path.exists(resultDir):
            os.mkdir(resultDir)

        #Spearman RMSD data
        fileNameList = ["spearmanRMSD.tsv", "Clustering.tsv", "spearmanCorrelations.tsv", "ABSpearmanCoeff.tsv", "knnOverlap.tsv", "PercentageCorrelationDetection.tsv"]
        for fileName in fileNameList:
            self.combine_files(newSimulationDir, resultDir, fileName, dublicationNum)
        

    def deleteExcessData(self, newSimulationDir, fileBenchmarksToKeep):
        folderName = os.path.basename(removesuffix(self.parameters.iniFile, '.ini'))
        dirName = newSimulationDir + "/" + folderName
        #if the file is not one of those three to keep, delete it
        if(not self.keepData and not( folderName.endswith(str(fileBenchmarksToKeep[0])) or folderName.endswith(str(fileBenchmarksToKeep[1])) or folderName.endswith(str(fileBenchmarksToKeep[2])) )):
            shutil.rmtree(dirName)

