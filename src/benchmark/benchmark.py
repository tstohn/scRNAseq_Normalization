#!/usr/bin/python3
from pickletools import string1
import sys
from datasetParsingFunctions import load_datasets_for_benchmark
from NormalizedDataHandler import NormalizedDataHandler
import argparse

def parse_args(argv):

    parser = argparse.ArgumentParser(description='Benchmarking all Normalization methods in a directory.')
    
    parser.add_argument('--groundtruth', help='compare norm methods to groundtruth. (Default simply benchmark norm methods: correlations, etc.)', action='store_true')
    parser.add_argument('--deleteOldData', help='instead of adding new data to the already existing folder in ./bin/BENCHMARK, the folder is deleted and a totally new Benchmark is created', action='store_true')
    
    parser.add_argument('--filterAbTypeForSpearmanCorrelation', help='run spearman correlation additionally for a subset of the data, when looking only at specific protein types',
                        type=str)
    parser.add_argument('dir', metavar='DIR', type=str)

    args = parser.parse_args()
    return(args)

def make_benchmark(dataset, groundtruth, deleteBenchmark, spearmanFilter):
    #initialization
    benchmark = NormalizedDataHandler(dataset, deleteBenchmark)

    #additional visualizations
    #benchmark.draw_tsne()
    #benchmark.calculate_MMD_between_treatments()

    #elaborate analysis of norm methods
    #benchmark.run_treatment_classification()
    benchmark.ab_spearman_correlation(groundtruth) # make RMSD of correlation differences => barplot

    if(groundtruth):
        #RMSD of real to norm data => barplot
        benchmark.validate_normalizedData_against_groundTruth()
        #calculate detected correlations of proteins
        #calculate removed batch effect
        #calculate detected treatment effect

    #additional correlation analysis for a subset od the data
    if(spearmanFilter):
        benchmark.ab_spearman_correlation(groundtruth, spearmanFilter)

def main():
    if(len(sys.argv) < 2):
        print("ERROR: use script with <python3 benchmark.py [directory of several datasets || directory with normalized files of one dataset]>\n")
        exit(1)
    print("Running benchmark of normalized scRNAseq data from: "+sys.argv[1]+"\n")

    args = parse_args(sys.argv)

    datasets = load_datasets_for_benchmark(args.dir)
    print("Found " + str(len(datasets)) + " datasets to run normalization benchmark on")

    #make classifications
    for dataset in datasets:
        make_benchmark(dataset, args.groundtruth, args.deleteOldData, args.filterAbTypeForSpearmanCorrelation)

if __name__ == '__main__':
    main()