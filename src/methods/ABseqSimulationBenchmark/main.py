import numpy as np
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

def parse_args():

    parser = argparse.ArgumentParser(description='Benchmarking Simulations')
    parser.add_argument('--t', help='threads',
                        type=int, default=1)
    parser.add_argument('dir', metavar='DIR', type=str)

    args = parser.parse_args()
    return(args)

def runSimulation(ini, newSimulationDir):
    try:
        param = Parameters(ini)
        benchmark = Benchmark(param)
        benchmark.run()
        benchmark.moveIntoOneFolder(newSimulationDir)
    except:
        print("Could not run Benchmark on " + ini)

def main():

    if(len(sys.argv) < 2):
        print("ERROR: use script with \'python3 main.py <folder_with_simulation.inis>\'\n")
        exit(1)

    args = parse_args()
    parameterFolder = args.dir

    #list all ini files
    iniFileList = [f for f in listdir(parameterFolder) if isfile(join(parameterFolder, f))]

    regex = re.compile('.*\.ini')
    iniFileList = [join(parameterFolder,x) for x in iniFileList if regex.match(x)]

    #generate a master directory for all experiments
    newSimulationDir = "./bin/BENCHMARKED_DATASETS/Simulations_" + str(datetime.now())
    os.makedirs(newSimulationDir)

    #run this as threads
    pool_size = args.t
    pool = Pool(pool_size)

    for ini in iniFileList:
        pool.apply_async(runSimulation, args=(ini, newSimulationDir))
            
    pool.close()
    pool.join()

if __name__ == '__main__':
    main()