# How To: single-cell protein simulations

#RUNNING A SIMULATION BENCHMARK:

look in Makefile for exmaples:
- e.g., run src/simulation/ABseqSimulationBenchmark/main.py <simulation.ini>

Ini examples are in simulation_ini.

The script will create three folders:
    bin/SIMULATIONS: this contains the raw simulated data
    bin/NORMALIZED_DATASETS: this contains the data normalized with the specified methods in the simulation.ini
    bin/BENCHMARKED_DATASETS: this will contain benchmarks on the various norm. methods, from correlation, to differential abundance, etc. ...

## Main parameters for simualtions:

*ProteinLevels*: a list of pairs [[A,B], [A,B], ...], where A is again a pair of μ, size which defines the neg.binomial distorbution
of protein μs. B is the number of proteins that should have a μ from this distribution. The final single-cell protein counts are then drawn from a neg.binomial dist with a μ from the above mention dist. and a size which is calculated with a function given by 'size' parameter to account for overdispersion in the data.
e.g.: [[[1500,2.0],50],...]

*size*: a function defining how the size parameter for the neg.binom. distributions for every protein count depends on μ to account for in overdispersion.

*libSize*: a pair defining mean and std. for log-normal distributions of sequencing depth effects. Every protein count for a cell will be DIVIDED by the corresponding lib-size factor drawn for that cell from this log-normal dist.

*noise*: random factor to multiply counts with, drawn with a gaussian dist. This parameters defines the std. of this gaussian dist.

BE AWARE: It is not possible to 'handpick' correlations in a correltion matrix, as this matrix must be positive semi-definite!!!
Therefore below parameters can be sampled, then we create the CLOSEST positive semi-definite matrix to these parameters. 
*proteinCorrelationMean* + *proteinCorrelationStd*: The mean and std. for a guassian dist. describing the average (background) protein correlations. For every protein pair we draw from this distribution.
*proteinCorrelationDists*: an array of triplets defining additional 'clusters' of protein correlations that we want to see in the data. E.g., [0.5, 0.1, 10] defines additional 10 proteins that should be inter-correlated by on average a correlation around 0.5 with a std. of 0.1.

If the above definition of correlationsdoes not work, you can use gaussian-copula to introduce correlations.
therefore you need to set *randomProteinCorrelations=1*. The amount of correlations in the data can then be controlled with a beta parameter: *betaParameter=x*
The beta parameter defines the amount of correlations (positive as well as negative).
smaller beta: stronger pos. and also negative correlations (beta must be > 0) 
higher beta: overall weaker correlations.
If you want to play around with this parameter as a starting point we suggest following: For 60 proteins a beta between 0.2 and 10 seems to generate very strong to very weak correlation matrices.
