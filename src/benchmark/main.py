#!/usr/bin/python3
from pickletools import string1
import sys

from numpy.lib.arraysetops import unique
from datasetParsingFunctions import load_datasets_for_benchmark
import argparse
import os 
import pandas as pd

sys.path.append('./src/simulation/ABseqSimulation')
sys.path.append('./src/methods/ToolBox')
from functions import *

'''
All datasets will be normalized, if a file with the same name already exists, it will be deleted before.
'''

def parse_args():

    parser = argparse.ArgumentParser(description='Benchmarking all Normalization methods in a directory.')
    
    parser.add_argument('--groundtruth', help='compare norm methods to groundtruth. (Default simply benchmark norm methods: correlations, etc.)', action='store_true')
    parser.add_argument('--deleteOldData', help='instead of adding new data to the already existing folder in ./bin/BENCHMARK, the folder is deleted and a totally new Benchmark is created', action='store_true')
    
    parser.add_argument('--iniFile', help='ini file - e.g. the same as from simulation, which states which proteins are correlated',
                        type=str)

    parser.add_argument('--filterAbTypeForSpearmanCorrelation', help='run spearman correlation additionally for a subset of the data, when looking only at specific protein types',
                        type=str)
    parser.add_argument('--stdout', help='write unimportant messages to a file', default="",
                        type=str)
    #thread is the number of parallel runs of benchmarks: however it does not work on Linux server, therefore if threads is -1, we run one benchmark after the other
    #but allow sklearn for treatment detection to spawn threads as it wants
    parser.add_argument('--t', help='threads',
                        type=int, default=-1)
    parser.add_argument('dir', metavar='DIR', type=str)

    args = parser.parse_args()
    return(args)

def make_benchmark(dataset, groundtruth, deleteBenchmark, spearmanFilter, iniFile, threadsSklearn):
    #initialization
    from NormalizedDataHandler import NormalizedDataHandler
    import Simulation

    benchmark = NormalizedDataHandler(dataset, groundtruth, deleteBenchmark, threadsSklearn)

    #additional visualizations
    benchmark.draw_tsne()
    benchmark.calculate_MMD_between_treatments()

    #elaborate analysis of norm methods
    treatmentNum = 1
    dataTmp = pd.read_csv(dataset[0], sep = '\t')
    treatmentNum = len(unique(dataTmp["cluster_id"]))
    
    if(treatmentNum > 1):
        try:
            benchmark.run_treatment_classification()
        except Exception as e: 
            print(e)
            printToTerminalOnce("\n ERROR: Treatment classification failed\n")

    if(groundtruth):
        params = Simulation.Parameters(iniFile)

        try:
            benchmark.ab_spearman_correlation(groundtruth) # make RMSD of correlation differences => barplot
        except Exception as e: 
            print(e)
            printToTerminalOnce("\n ERROR: Spearman correlation failed\n")
        #RMSD of fold cahnges between total AB counts per sample 
        #(idea: insample fold cahnges between dofferent protein counts stay the same after normalization, only the different samples are scaled)
        try:
            benchmark.validate_normalizedData_against_groundTruth()
        except Exception as e: 
            print(e)
            printToTerminalOnce("\n ERROR: RMSD between ABcount to min ABcount failed\n")

        #calculate detected correlations of proteins - check we have wanted and not unwanted corr
        #make heatmap of all wanted, and of 50 randomly chosen unwanted variations
        try:
            benchmark.validate_correlations(params.proteinCorrelations)
        except Exception as e: 
            print(e)
            printToTerminalOnce("\n ERROR: Detection of wanted Correlation failed\n")

        if(params.diffExProteins != None):
            try:
                benchmark.validate_treatment_effect(params.diffExProteins, params.treatmentVector)
            except Exception as e: 
                print(e)
                printToTerminalOnce("\n ERROR: Treatment Effect validation failed\n")

        #calculate detected treatment effect - check we have wanted treatment effect and not unwanted
        if(params.batchFactors != None):
            #calculate removed batch effect
            try:
                benchmark.validate_batch_effect()
            except Exception as e: 
                print(e)
                printToTerminalOnce("\n ERROR: Batch Eff failed\n")

        #additional correlation analysis for a subset od the data
        if(spearmanFilter):
            try:
                benchmark.ab_spearman_correlation(groundtruth, spearmanFilter)
            except Exception as e: 
                print(e)
                printToTerminalOnce("\n ERROR: Spearman Correlelation for proteins of type " + spearmanFilter + " failed\n")
    else:

        #SPEARMAN CORRELATION WITHOUT GROUNDTRUTH
        try:
            benchmark.ab_spearman_correlation() # make RMSD of correlation differences => barplot
        except Exception as e: 
            print(e)
            printToTerminalOnce("\n ERROR: Spearman correlation failed\n")

def main():
    if(len(sys.argv) < 2):
        print("ERROR: use script with <python3 benchmark.py [directory of several datasets || directory with normalized files of one dataset]>\n")
        exit(1)

    args = parse_args()
    if(args.stdout != ""):
        outfile = open(args.stdout, 'a+')
        sys.stdout = outfile
        sys.stderr = outfile
    threadsSklearn = args.t
    #necessay for pool to work with number of threads on linux: Issue: https://github.com/numpy/numpy/issues/14474
    
    if(threadsSklearn != -1):
        os.environ['OPENBLAS_NUM_THREADS'] = str(threadsSklearn)
        os.environ["OMP_NUM_THREADS"] = str(threadsSklearn)
        os.environ["MKL_NUM_THREADS"] = str(threadsSklearn)
        os.environ["VECLIB_MAXIMUM_THREADS"] = str(threadsSklearn)
        os.environ["NUMEXPR_NUM_THREADS"] = str(threadsSklearn)
        threadsSklearn = 1
    printToTerminalOnce("Running benchmark of normalized scRNAseq data from: "+sys.argv[len(sys.argv)-1]+"\n")
    datasets = load_datasets_for_benchmark(args.dir)
    printToTerminalOnce("Found " + str(len(datasets)) + " datasets to run normalization benchmark on")

    #make classifications
    for dataset in datasets:
        make_benchmark(dataset, args.groundtruth, args.deleteOldData, args.filterAbTypeForSpearmanCorrelation, args.iniFile, threadsSklearn)

    if(args.stdout != ""):
        outfile.close()

if __name__ == '__main__':
    main()