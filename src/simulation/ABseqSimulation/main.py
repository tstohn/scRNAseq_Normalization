from Simulation import SingleCellSimulation, Parameters
import numpy as np
import sys
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Simulations of ABCount Data')
    parser.add_argument('--stdout', help='write unimportant messages to a file', default="",
                        type=str)
    parser.add_argument('dir', metavar='DIR', type=str)

    args = parser.parse_args()
    return(args)

def main():

    args = parse_args()
    parameter_file = args.dir
    if(args.stdout != ""):
        outfile = open(args.stdout, 'a+')
        sys.stdout = outfile
        sys.stderr = outfile

    print("Starting Simulation:")
    param = Parameters(parameter_file)
    sim = SingleCellSimulation(param)
    sim.simulateData()
    sim.save_data()
    outfile.close()

if __name__ == '__main__':
    main()