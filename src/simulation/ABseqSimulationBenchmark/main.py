import sys
from simBenchmark import Parameters,Benchmark
from os import listdir
from os.path import isfile, join
import regex as re
import os
from datetime import datetime  
import shutil
import argparse
from multiprocessing.pool import ThreadPool as Pool
sys.path.append('./src/methods/ToolBox')
from functions import *

"""
    Default Usage: IniFile with one line that contains the string INIRANGE: this line defines the parameter that is altered fir simulation purposes

    Parameters:
        inputFile: the input parameter file, it can contain ranges for certain values, for each value in this range a new ini file is temporaryly generated
             to store the parameters during simulation, since several different tools will run and need access to it.
             For each Simulation we can have ONLY ONE variable parameter, that is read from the file and for each variable value a new tmp file
             is generated.
             By default all normalized and simulated datasets are deleted after simulation to not overflow in data (the ini/ tsv file of simulated data
             is also copied into the folder for normalizations - those files are deleted as well)
        -s: add flag <-s> to use a directory instead of a single ini files with ranges
            if s variable is set we have a SINGLE simulation file for all simulations (default it is set), if false we have several inis and in this case -in MUST BE A DIRECTORY
        -k : if set we keep all simulation and normalization data, as well as simulated inis
        -b : to run benchmark only for a bunch of files (e.g. if benchmark failed but simulations/ normalizaiton ran perfectly).
             This means we MUST HAVE a DIRECTORY WITH ALL INIs and a several DIRECTORIES FOR EVERY INI WITH NORMALIZED FILES in bin
        --d : how often an experiemnt is repeated and averaged, on average this is 5 times
        --t : threads
        
        -c : subset the cell number, and run it as a new benchmark, this keeps exaclty the same cells and runs the benchmark
"""
def parse_args():

    parser = argparse.ArgumentParser(description='Benchmarking Simulations')
    #if minus one, do not explicitely restrict threads - allowing python, R to spwan threads in their libraries
    parser.add_argument('--t', help='threads',
                        type=int, default=-1)
    parser.add_argument('dir', metavar='DIR', type=str)
    parser.add_argument('-s', action='store_false')
    parser.add_argument('-b', action='store_true')
    parser.add_argument('-k', action='store_true')
    parser.add_argument('-c', action='store_true')
    parser.add_argument('--d', help='duplicates',
                        type=int, default=5)
    
    #variables for the benchmarking for classification/ knn Overlap (for now we always apply zscore)
    parser.add_argument('--knn', help='number of KNN to calculate overlap of groundtruth/ normalized neighbors & for classification (default 20, -1 to skip step, 0 for a gradient of K which is only used for overlap)', type=int, default=20)
    parser.add_argument('--metric', help = 'metric used for knn overlap/ classification: euclidean, manhattan or cosine',
                        type=str, default='manhattan')

    args = parser.parse_args()
    
    return(args)

def map_three_ini_values(info, list):
    elementNum = 0
    for element in info:
        if(elementNum == 0):
            list[0] = float(element)
        elif(elementNum == 1):
            list[1] = float(element)
        elif(elementNum == 2):
            list[2] = float(element)
        elementNum += 1

#the newly generated files are written into a tmp directory, there they consist of a name which is basically a number from 0 to <numberOfSimulatedSamples>
def generate_simulation_iniFiles(iniFile, fileBenchmarksToKeep):
    import numpy as np
    if(os.path.isdir(os.path.realpath(iniFile))):
        printToTerminalOnce("ERROR: Directiory given, but the flag for handling single predefined ini files not set; ABORT SIMULATION")
        exit()
    printToTerminalOnce("Generating INI Files for simulations")
    dir_path = os.path.dirname(os.path.realpath(iniFile)) + "/TmpIniDir/"
    if(os.path.exists(dir_path)):
        shutil.rmtree(dir_path)
    os.mkdir(dir_path)

    #reading different range parameters and creating files for those
    variableParameter = ""
    file = open(iniFile, "r")
    line = file.readline()
    start = end = factor = 1
    start2 = end2 = factor2 = 1
    value2Set = False
    foundIni = False

    #INIRANGE IS ALWAYS READ THE SAME WAY:
    # we have three comma seperated values that describe a range,
    # we might also have a second range that is seperated by a semicolon
    while line:
        if( (not line.startswith("#")) and ("INIRANGE" in line) ):
            info = re.match(("(.+?)INIRANGE=(.*)"), line)

            variableParameter = str(info[1]) 
            infoArray = str(info[2]).split(";")
            #if the range itself has several values (so far only 2) parse them seperately
            if(len(infoArray) > 1):
                valueList_1 = [start, end, factor]
                valueList_2 = [start, end, factor]

                info = str(infoArray[0]).split(",")
                map_three_ini_values(info, valueList_1)
                info = str(infoArray[1]).split(",")
                map_three_ini_values(info, valueList_2)

                start = valueList_1[0]
                end = valueList_1[1]
                factor = valueList_1[2]

                start2 = valueList_2[0]
                end2 = valueList_2[1]
                factor2 = valueList_2[2]
                value2Set = True
            else:
                info = str(info[2]).split(",")
                valueList = [start, end, factor]
                map_three_ini_values(info, valueList)
                start = valueList[0]
                end = valueList[1]
                factor = valueList[2]

            if(foundIni==True):
                printToTerminalOnce("ERROR parsing INIRANGE for Benchmark; foundn several lines with 'INIRANGE', but a simulation can only have a single variable\n")
                exit()
            foundIni = True
        line = file.readline()
    file.close()
    
    count = 0
    fileBenchmarksToKeep.append(0)

    #if we have no INIRANGE: write at least the file itself by setting end to 2
    #to initiate a single iteration
    if(start == 1 and end == 1 and factor == 1):
        end = 2
    for i in np.arange(start, end, factor):
        j = 0
        if value2Set:
            j = start2 + count*factor2

        if(len(fileBenchmarksToKeep)==1 and count>=((end-start)/factor)/2):
            fileBenchmarksToKeep.append(count)
        i = round(i,2)
        newFile = open(dir_path + "/" + str(count) + ".ini", "a")
        file = open(iniFile, "r")
        line = file.readline()
        while line:
            if(str.startswith(line, variableParameter+"INIRANGE") and (variableParameter != "")):
                line = file.readline()
                continue
            elif(str.startswith(line, variableParameter) and (variableParameter != "")):
                #write new variable line

                #LIBRARY SIZE
                if(variableParameter=="libSize"):
                    newLine = line.replace("X", str(i))
                    newFile.write(newLine)
                    print(newLine)

                #NOISE
                if(variableParameter=="noise"):
                    newFile.write("noise=" + str(i) + "\n")
                if(variableParameter=="ProteinLevels"):
                    integerProteinCount = int(i)
                    newLine = line.replace("X", str(integerProteinCount))
                    if value2Set:
                        integerProteinCount2 = int(j)
                        newLine = newLine.replace("Y", str(integerProteinCount2))
                    newFile.write(newLine)
                if(variableParameter=="proteinNoise"):
                    newFile.write("proteinNoise=" + str(i) + "\n")
                
                #CORRELATION LINE
                #same as ProteinLevels: we can have diff distributions of correlations
                #think of it as the correlation density plots: now we define how many protein-pairs should be sampled
                #from those distributions that we define
                if(variableParameter=="proteinCorrelationDists"):
                    proteinCorrParam = float(i)
                    newLine = line.replace("X", str(proteinCorrParam))
                    if value2Set:
                        proteinCorrParam2 = float(j)
                        newLine = newLine.replace("Y", str(proteinCorrParam2))
                    newFile.write(newLine)

                #CLUSTER LINE
                #increasing number of proteins effected by cluster (might make most sense to set numberProteins=[X] with ONE additionalCluster, but other scenarios are possible)
                if(variableParameter=="numberProteins"):
                    newLine = line.replace("X", str(i))
                    newFile.write(newLine)
                #increasing factor weight to scale prrteins between clsuters (might make most sense to set abundanceFactors=[X] with ONE additionalCluster, but other scenarios are possible)
                if(variableParameter=="abundanceFactors"):
                    newLine = line.replace("X", str(i))
                    newFile.write(newLine)
                #increasing percentage for correlationFactor (might make most sense to set correlationFactors=[X] with ONE additionalCluster, but other scenarios are possible)
                if(variableParameter=="correlationFactors"):
                    newLine = line.replace("X", str(i))
                    newFile.write(newLine)
                if(variableParameter=="betaParameter"):
                    newLine = line.replace("X", str(i))
                    newFile.write(newLine)  
                if(variableParameter=="betaFactors"):
                    newLine = line.replace("X", str(i))
                    newFile.write(newLine)  
                #increasing cellPercentages (have in mind, all NON-X cell percentages get evenly scaled to sum to 100)
                if(variableParameter=="cellPercentages"):
                    #we need to scale other cell population percentages so that they sum to 100
                    info = re.match(("cellPercentages=\[(.*)\]"), line)
                    info = str(info[1]).split(",")
                    fillTo = 100.0 - float(i)
                    recentValue = 0.0
                    difference = 0.0
                    cellPercentagesLength = len(info) - 1
                    for el in info:
                        if(el != "X"):
                            recentValue += float(el)
                    if(recentValue is not fillTo):
                        difference = (fillTo - recentValue)
                    newLine = "cellPercentages=["
                    differenceReached = 0
                    infoIdx = 0
                    cellPercentageList = [element for element in info if element != "X"]
                    for el in cellPercentageList:
                        newValue = 0.0
                        if(el == "X"):
                            continue
                        if(infoIdx == cellPercentagesLength-1):
                            newValue = float(el) + (difference - differenceReached)
                        else:
                            newValue = float(el) + round((difference/cellPercentagesLength),1)
                            differenceReached += round((difference/cellPercentagesLength),1)
                        infoIdx += 1
                        newLine += str(newValue) + ","
                    newLine += str(i) + "]\n"
                    newFile.write(newLine)
                if(variableParameter=="numberClusterSpecificProteins"):
                    integerProteinCount = int(i)
                    newLine = line.replace("X", str(integerProteinCount))
                    newFile.write(newLine)
            else:
                #write same line into file
                newFile.write(line)
            line = file.readline()
        count += 1
        newFile.close()
        file.close()
    fileBenchmarksToKeep.append(count-1) #keep the number of the last file
    return(dir_path)

def delete_tmp_folder(folder):
    if(os.path.exists(folder)):
        shutil.rmtree(folder)

def runSimulation(ini, benchmarkBaseDir, stdoutFile, noExplicitelySetThreads, duplicates, keepData, fileBenchmarksToKeep, 
                  recentSimulationNumber, numberOfSimulations, benchmarkOnly, knnOverlap, knnMetric, subsetCells):
    printToTerminalOnce("#  RUNNING SIMULATION[" + str(recentSimulationNumber) + "/" + str(numberOfSimulations) + "]")
    try:
        for duplicateIdx in range(0,duplicates):
            param = Parameters(ini)
            benchmark = Benchmark(param, benchmarkBaseDir, stdoutFile, noExplicitelySetThreads, keepData, benchmarkOnly, knnOverlap, knnMetric, subsetCells, duplicateIdx)
            benchmark.run()
            
            #copy all the individual results in bin/BENCHMARK/xyz into one shared directory 
            #and create one summary file
            benchmark.copyResultsIntoOneFile(benchmarkBaseDir, duplicateIdx)
            #if not keepDta==TRUE we only keep three representative datasets for benchmarks (first, middle and last one)
            #it also only keeps the datasets of all cells (not for subsets of cells)
            benchmark.deleteExcessData(benchmarkBaseDir, fileBenchmarksToKeep) 
            
    except Exception as e: 
        print("#ERROR MESSAGE: \'" + str(e) + "\'\n")
        printToTerminalOnce("\n ERROR: Could not run Benchmark on " + ini + "\n")

def main():

    if(len(sys.argv) < 2):
        print("ERROR: use script with \'python3 main.py <iniFile.tsv>\'\n")
        exit(1)

    stdoutFile = 'LogFile.txt'
    outfile = open(stdoutFile, 'w')
    sys.stdout = outfile

    sys.stderr = outfile
    args = parse_args()

    try:
        #set thread limit before importing numpy
        pool_size = args.t
        noExplicitelySetThreads = False
        if(pool_size == -1):
            pool_size = 1
            noExplicitelySetThreads = True
        else:
            #necessay for pool to work with number of threads on linux: Issue: https://github.com/numpy/numpy/issues/14474
            os.environ['OPENBLAS_NUM_THREADS'] = str(pool_size)
            os.environ["OMP_NUM_THREADS"] = str(pool_size)
        import numpy as np


        tmpDir = ""
        keepData = False
        benchmarkOnly = False
        subsetCells = False

        fileBenchmarksToKeep = []
        if(args.s):
            tmpDir = generate_simulation_iniFiles(args.dir, fileBenchmarksToKeep)
        else:
            tmpDir = args.dir
        if(args.k):
            keepData = True
        if(args.b):
            benchmarkOnly = True
        if(args.c):
            subsetCells = True

        parameterFolder = tmpDir

        #list all ini files
        iniFileList = [f for f in listdir(parameterFolder) if isfile(join(parameterFolder, f))]
        
        regex = re.compile('.*\.ini')
        iniFilePathList = [join(str(parameterFolder),x) for x in iniFileList if regex.match(x)]
        #generate a master directory for all experiments
        simulationName = os.path.splitext(os.path.basename(args.dir))[0]
        benchmarkBaseDir = "./bin/BENCHMARKED_DATASETS/Simulations_" + simulationName + "_" + str(datetime.now())
        os.makedirs(benchmarkBaseDir)

        #run this as threads
        pool = Pool(pool_size)
        numberOfSimulations = len(iniFilePathList)
        recentSimulationNumber = 1
        for ini in iniFilePathList:
            pool.apply_async(runSimulation, args=(ini, benchmarkBaseDir, stdoutFile, noExplicitelySetThreads, args.d, keepData, fileBenchmarksToKeep, recentSimulationNumber, numberOfSimulations, benchmarkOnly, args.knn, args.metric, subsetCells))
            recentSimulationNumber = recentSimulationNumber + 1

        pool.close()
        pool.join()
        outfile.close()

    except Exception as e: 
        printToTerminalOnce("Benchmark Failed")
        print("#ERROR MESSAGE: \'" + str(e) + "\'\n")

    if(args.s and not args.k):
        delete_tmp_folder(tmpDir)

if __name__ == '__main__':
    main()