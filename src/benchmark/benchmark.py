#!/usr/bin/python3
import sys
from datasetParsingFunctions import load_datasets_for_benchmark
from NormalizedDataHandler import NormalizedDataHandler

def make_benchmark(dataset):
    benchmark = NormalizedDataHandler(dataset)
    benchmark.dt_classification()
    benchmark.knn_clasification()
    benchmark.draw_tsne()
    benchmark.ab_spearman_correlation()

    #dataset scpecific analyses
    if(benchmark.dataset_name.startswith("scIDseq-vanEijl")):
        benchmark.ab_spearman_correlation("total")

def main():
    if(len(sys.argv) != 2):
        print("ERROR: use script with <python3 benchmark.py [directory of several datasets || directory with normalized files of one dataset]>\n")
        exit(1)
    print("Running benchmark of normalized scRNAseq data from: "+sys.argv[1]+"\n")

    datasets = load_datasets_for_benchmark(sys.argv[1])
    print("Found " + str(len(datasets)) + " datasets to run normalization benchmark on")

    #make classifications
    for dataset in datasets:
        make_benchmark(dataset)

if __name__ == '__main__':
    main()