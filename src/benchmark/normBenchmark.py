#!/usr/bin/python3
from pickletools import string1
import sys
from datasetParsingFunctions import load_datasets_for_benchmark
from NormalizedDataHandler import NormalizedDataHandler
import argparse

sys.path.append('./src/methods/ABseqSimulation')
import Simulation

def parse_args(argv):

    parser = argparse.ArgumentParser(description='Benchmarking all Normalization methods in a directory.')
    
    parser.add_argument('--groundtruth', help='compare norm methods to groundtruth. (Default simply benchmark norm methods: correlations, etc.)', action='store_true')
    parser.add_argument('--deleteOldData', help='instead of adding new data to the already existing folder in ./bin/BENCHMARK, the folder is deleted and a totally new Benchmark is created', action='store_true')
    
    parser.add_argument('--iniFile', help='ini file - e.g. the same as from simulation, which states which proteins are correlated',
                        type=str)

    parser.add_argument('--filterAbTypeForSpearmanCorrelation', help='run spearman correlation additionally for a subset of the data, when looking only at specific protein types',
                        type=str)
    parser.add_argument('dir', metavar='DIR', type=str)

    args = parser.parse_args()
    return(args)

def make_benchmark(dataset, groundtruth, deleteBenchmark, spearmanFilter, iniFile):
    #initialization
    benchmark = NormalizedDataHandler(dataset, groundtruth, deleteBenchmark)

    #additional visualizations
    benchmark.draw_tsne()
    benchmark.calculate_MMD_between_treatments()

    #elaborate analysis of norm methods
    try:
        benchmark.run_treatment_classification()
    except:
       print("Treatment classification failed")
    try:
        benchmark.ab_spearman_correlation(groundtruth) # make RMSD of correlation differences => barplot
    except:
        print("Spwearman correlation failed")

    if(groundtruth):
        params = Simulation.Parameters(iniFile)

        #RMSD of fold cahnges between total AB counts per sample 
        #(idea: insample fold cahnges between dofferent protein counts stay the same after normalization, only the different samples are scaled)
        try:
            benchmark.validate_normalizedData_against_groundTruth()
        except:
            print("RMSD between ABcount to min ABcount failed")


        #calculate detected correlations of proteins - check we have wanted and not unwanted corr
        #make heatmap of all wanted, and of 50 randomly chosen unwanted variations
        try:
            benchmark.validate_correlations(params.proteinCorrelations)
        except:
            print("Detection of wanted Correlation failed")

    #calculate detected treatment effect - check we have wanted treatment effect and not unwanted
    if(params.batchFactors != None):
        #calculate removed batch effect
        try:
            benchmark.validate_batch_effect()
        except:
            print("Batch Eff failed")
    #additional correlation analysis for a subset od the data
    if(spearmanFilter):
        try:
            benchmark.ab_spearman_correlation(groundtruth, spearmanFilter)
        except:
            print("Spearman Correlelation for proteins of type " + spearmanFilter + " failed")
    if(params.diffExProteins != None):
        try:
            benchmark.validate_treatment_effect(params.diffExProteins, params.treatmentVector)
        except:
            print("Treatment Effect validation failed")

def main():
    if(len(sys.argv) < 2):
        print("ERROR: use script with <python3 benchmark.py [directory of several datasets || directory with normalized files of one dataset]>\n")
        exit(1)
    print("Running benchmark of normalized scRNAseq data from: "+sys.argv[len(sys.argv)-1]+"\n")

    args = parse_args(sys.argv)

    datasets = load_datasets_for_benchmark(args.dir)
    print("Found " + str(len(datasets)) + " datasets to run normalization benchmark on")

    #make classifications
    for dataset in datasets:
        make_benchmark(dataset, args.groundtruth, args.deleteOldData, args.filterAbTypeForSpearmanCorrelation, args.iniFile)

if __name__ == '__main__':
    main()