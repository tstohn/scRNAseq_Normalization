from Simulation import SingleCellSimulation, Parameters
import numpy as np
import sys

def main():

    if(len(sys.argv) != 2):
        print("ERROR: use script with \'python3 main.py <simulation_settings.ini>\'\n")
        exit(1)
    parameter_file = sys.argv[1]

    param = Parameters(parameter_file)
    sim = SingleCellSimulation(param)
    sim.simulateData()
    #sim.save_data(groundTruth=True)
    #sim.save_data()

if __name__ == '__main__':
    main()