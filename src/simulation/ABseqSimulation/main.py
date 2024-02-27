import sys
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Simulations of ABCount Data')
    parser.add_argument('--stdout', help='write unimportant messages to a file', default="",
                        type=str)
    parser.add_argument('dir', metavar='DIR', type=str)
    parser.add_argument('--t', help='threads', type=int, default=-1)

    args = parser.parse_args()
    return(args)

def main():
    print("Starting Simulation:\n")

    args = parse_args()
    parameter_file = args.dir
    if(args.stdout != ""):
        outfile = open(args.stdout, 'a+')
        sys.stdout = outfile
        sys.stderr = outfile

    #set thread environment variables and load libraries
    import os
    if(args.t != -1):
        os.environ["OMP_NUM_THREADS"] = str(args.t) # export OMP_NUM_THREADS=4
        os.environ["OPENBLAS_NUM_THREADS"] = str(args.t) # export OPENBLAS_NUM_THREADS=4 
        os.environ["MKL_NUM_THREADS"] = str(args.t) # export MKL_NUM_THREADS=6
        os.environ["VECLIB_MAXIMUM_THREADS"] = str(args.t) # export VECLIB_MAXIMUM_THREADS=4
        os.environ["NUMEXPR_NUM_THREADS"] = str(args.t) # export NUMEXPR_NUM_THREADS=6
    from Simulation import SingleCellSimulation, Parameters #library imports numpy, need to set environment variables first

    param = Parameters(parameter_file)
    sim = SingleCellSimulation(param)
    sim.simulateData()
    sim.save_data()

    if(args.stdout != ""):
        outfile.close()

if __name__ == '__main__':
    main()