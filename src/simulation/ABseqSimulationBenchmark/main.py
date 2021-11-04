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

def parse_args():

    parser = argparse.ArgumentParser(description='Benchmarking Simulations')
    parser.add_argument('--t', help='threads',
                        type=int, default=1)
    parser.add_argument('dir', metavar='DIR', type=str)

    args = parser.parse_args()
    return(args)

def runSimulation(ini, newSimulationDir, stdoutFile):
    try:
        param = Parameters(ini)
        benchmark = Benchmark(param, stdoutFile)
        benchmark.run()
        benchmark.moveIntoOneFolder(newSimulationDir)
        benchmark.moveIntoOneFolder(stdoutFile)
    except:
        printToTerminalOnce("\n ERROR: Could not run Benchmark on " + ini + "\n")

def main():

    if(len(sys.argv) < 2):
        print("ERROR: use script with \'python3 main.py <folder_with_simulation.inis>\'\n")
        exit(1)
    stdoutFile = 'LogFile.txt'
    outfile = open(stdoutFile, 'w')
    sys.stdout = outfile
    sys.stderr = outfile

    args = parse_args()
    #set thread limit before importing numpy
    pool_size = args.t
    #necessay for pool to work with number of threads on linux: Issue: https://github.com/numpy/numpy/issues/14474
    os.environ['OPENBLAS_NUM_THREADS'] = str(pool_size)
    import numpy as np

    parameterFolder = args.dir

    #list all ini files
    iniFileList = [f for f in listdir(parameterFolder) if isfile(join(parameterFolder, f))]

    regex = re.compile('.*\.ini')
    iniFileList = [join(parameterFolder,x) for x in iniFileList if regex.match(x)]

    #generate a master directory for all experiments
    newSimulationDir = "./bin/BENCHMARKED_DATASETS/Simulations_" + str(datetime.now())
    os.makedirs(newSimulationDir)

    #run this as threads
    pool = Pool(pool_size)
    for ini in iniFileList:
        pool.apply_async(runSimulation, args=(ini, newSimulationDir, stdoutFile))
            
    pool.close()
    pool.join()
    outfile.close()
    
if __name__ == '__main__':
    main()