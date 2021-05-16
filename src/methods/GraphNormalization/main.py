from NormalizationGraph import NormalizationGraph
import sys

def normalize_data(dataset):
    NormalizationGraph(dataset)

def main():
    if(len(sys.argv) != 2):
        print("ERROR: use script with <python3 graphNormalization.py [directory of datasets]>\n")
        exit(1)

    normalize_data(sys.argv[1])

if __name__ == '__main__':
    main()
