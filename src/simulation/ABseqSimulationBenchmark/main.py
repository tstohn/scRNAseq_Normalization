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
        -in: the input parameter file, it can contain ranges for certain values, for each value in this range a new ini file is temporaryly generated
             to store the parameters during simulation, since several different tools will run and need access to it.
             For each Simulation we can have ONLY ONE variable parameter, that is read from the file and for each variable value a new tmp file
             is generated.
             By default all normalized and simulated datasets are deleted after simulation to not overflow in data (the ini/ tsv file of simulated data
             is also copied into the folder for normalizations - those files are deleted as well)
        -s: if set we have a SINGLE simulation for each file in the directory, in this case -in MUST BE A DIRECTORY
        -k : if set we keep all simulation and normalization data
"""
def parse_args():

    parser = argparse.ArgumentParser(description='Benchmarking Simulations')
    #if minus one, do not explicitely restrict threads - allowing python, R to spwan threads in their libraries
    parser.add_argument('--t', help='threads',
                        type=int, default=-1)
    parser.add_argument('dir', metavar='DIR', type=str)
    parser.add_argument('-s', action='store_false')
    parser.add_argument('-k', action='store_true')


    args = parser.parse_args()
    
    return(args)

def generate_simulation_iniFiles(iniFile):
    print("generating")
    if(os.path.isdir(os.path.realpath(iniFile))):
        printToTerminalOnce("ERROR: Directiory given, but the flag for handling single predefined ini files not set; ABORT SIMULATION")
        exit()
    print("Generating INI Files for simulations")
    dir_path = os.path.dirname(os.path.realpath(iniFile)) + "/TmpIniDir/"
    if(os.path.exists(dir_path)):
        shutil.rmtree(dir_path)
    os.mkdir(dir_path)

    #reading different range parameters and creating files for those
    variableParameter = ""
    file = open(iniFile, "r")
    line = file.readline()
    start = end = factor = 1
    while line:
        line = file.readline()
        if( (not line.startswith("#")) and ("INIRANGE" in line) ):
            info = re.match(("(.+?)INIRANGE=(.*)"), line)
            variableParameter = str(info[1])        
            info = str(info[2]).split(",")
            elementNum = 0
            for element in info:
                if(elementNum == 0):
                    start = int(element)
                elif(elementNum == 1):
                    end = int(element)
                elif(elementNum == 2):
                    factor = int(element)
                elementNum += 1
            break
    file.close()
    
    count = 0
    for i in range(start, end, factor):
        newFile = open(dir_path + "/" + str(count) + ".ini", "a")
        file = open(iniFile, "r")
        line = file.readline()
        while line:
            line = file.readline()
            if(str.startswith(line, variableParameter+"INIRANGE")):
                continue
            elif(str.startswith(line, variableParameter)):
                #write new variable line
                if(variableParameter=="libSize"):
                    newFile.write("libSize=[1," + str(i) + "]\n")
            else:
                #write same line into file
                newFile.write(line)
        count += 1
    return(dir_path)

def delete_tmp_folder(folder):
    if(os.path.exists(folder)):
        shutil.rmtree(folder)

def runSimulation(ini, newSimulationDir, stdoutFile, noExplicitelySetThreads, keepData):
    try:
        param = Parameters(ini)
        benchmark = Benchmark(param, stdoutFile, noExplicitelySetThreads, keepData)
        benchmark.run()
        benchmark.moveIntoOneFolder(newSimulationDir)
        benchmark.moveIntoOneFolder(stdoutFile)
    except Exception as e: 
        print(e)
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

    tmpDir = ""
    keepData = False
    if(args.s):
        tmpDir = generate_simulation_iniFiles(args.dir)
    else:
        tmpDir = args.dir
    if(args.k):
        keepData = True

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

        parameterFolder = tmpDir

        #list all ini files
        iniFileList = [f for f in listdir(parameterFolder) if isfile(join(parameterFolder, f))]

        regex = re.compile('.*\.ini')
        iniFilePathList = [join(str(parameterFolder),x) for x in iniFileList if regex.match(x)]
        #generate a master directory for all experiments
        newSimulationDir = "./bin/BENCHMARKED_DATASETS/Simulations_" + str(datetime.now())
        os.makedirs(newSimulationDir)

        #run this as threads
        pool = Pool(pool_size)
        for ini in iniFilePathList:
            pool.apply_async(runSimulation, args=(ini, newSimulationDir, stdoutFile, noExplicitelySetThreads, keepData))
                
        pool.close()
        pool.join()
        outfile.close()
    except:
        printToTerminalOnce("Benchmark Failed")

    if(args.s):
        delete_tmp_folder(tmpDir)
    
if __name__ == '__main__':
    main()