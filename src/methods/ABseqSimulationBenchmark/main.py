import numpy as np
import sys
from Benchmark import Parameters,Benchmark
from os import listdir
from os.path import isfile, join
import regex as re
import os
from datetime import datetime  
import shutil

def main():

    if(len(sys.argv) != 2):
        print("ERROR: use script with \'python3 main.py <folder_with_simulation.inis>\'\n")
        exit(1)
    parameterFolder = sys.argv[1]

    #list all ini files
    iniFileList = [f for f in listdir(parameterFolder) if isfile(join(parameterFolder, f))]

    regex = re.compile('.*\.ini')
    iniFileList = [join(parameterFolder,x) for x in iniFileList if regex.match(x)]

    #generate a master directory for all experiments
    newSimulationDir = "./bin/BENCHMARKED_DATASETS/Simulations_" + str(datetime.now())
    os.makedirs(newSimulationDir)

    #run this as threads
    for ini in iniFileList:
        param = Parameters(ini)
        benchmark = Benchmark(param)
        benchmark.run()
        benchmark.moveIntoOneFolder(newSimulationDir)

if __name__ == '__main__':
    main()