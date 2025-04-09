# ScNormalisation


## TO DO:
- local venv in MAC is corrupted
CHANGE CORRELATION / SAMPLING FROM DISTRIBUTIONS

## DISTRIBUTIONS OF ProteinCounts/ Library-Size
- Splatter: mean of all genes from gamma dist, then poisson for cells variability:
    - for now stay with Neg.Binom.: maybe argue somewhere why we used it
- understand gaussian copula
- Libsize: from log-normal

## DATASETS
- REVII44_ITGB_countstbl_allinfo.csv: origional scIDseq datasets from paper (ITGB1 selected cells)
- REVII44_ag1478_countstbl.csv: scIDseq dataset of EGFR inh treated cells
- full_object.rds: phoshoseq data, containing 64 ABs from 4 donors (Fig2)

## REQUIREMENTS
make file contains two commands: 
update_requirements and install_requirements to handle python (venv) and R requirements (R libs not in anenvironment - script simply collects all lib names from scripts without version)
- make requirements file for python ( python3 -m pip freeze)
- PYTHON:         python3 -m pip install -r requirements
- R requirements have to still be installed manually

## INI SET-UP
**Settigng Up Pipeline**:
- settings.ini in git directory (write dorectory for files (sc*proteins) that should be normiazed - for simple normalize runs without whole simulation)
- this folder stores all raw-files that should be normalized, 

**Normalization INI**
simply maps column names of file to column-names as they r used in the scripts (normally not necessary if files r processed so that column names match)


**Simulation_ini file**:
an exmaple is in INI_EXAMPLES, with definitions of all veriables that can be set

## VARIABLES (INIRANGEs)
Library Size effect:
- only libSize parameter is used, but theoretically you can enable also pcrCycles, pcrCapture, seqAmplificationEfficiency: but they r commented out for now...

INIRANGE: the varibale parameter must have  two values X,Y that can be set for any arbitrary parameter from the 3 INIRANGE values
first value is X, second is Y, so you can use them to e.g. select different protein-ranges to sample, or use simply X to set different means for neg.binom distribution
INIRANGE: it consists of 3 parameters, start, end, steps (where end is exclusive: while(start + steps < end))
INIRANGE works with following parameter: libSize, noise, proteinNoise, proteinCorrelationFactors, proteinCorrelation,diffExProteins

#CORRELATIONS 
#INIRANGE:
#1.) we can simulate ranges of mean/std for correlations, e.g., increasing correlations between proteinss
defines sets of correlated proteins, where every set is defined by a gaussian distrinution from which correlations are drawn, sets are all INTERCONNECTED,
#proteins is really the number of proteins that r interconnected!! So with correlations we can ONLY SAMPLE SETS OF CORRELATED PROTEINS, it s in a way
not ONLY the number of correlations, but number of correlated sets. The sets select proteins randomly that r included, however, proteins DO NOT OVERLAP between sets - sets are UNIQUE sets of proteins.
So we have length(proteinCorrelationDists) + 1 distributions for correltions to sample from (the baseline dist around zero plus all additionally mentioned dist here).
We can have as many proteinCorrelationDists as we want.
#proteinCorrelationDists: \[[mean, sd, #proteins];[...]]
#proteinCorrelationMean/proteinCorrelationStd for all the other proteins: their correlations are drawns from this baseline distributions (default: mean=0, sstd=0.4)
#roteinCorrelationsINIRANGE: can be used to change mean, STD over a range, or two different means of two distr. of correlations, etc. 

#Cluster
When we simualte several clusters we basically simulate sets of proteins that have an increased/ decreased proteins count. Additioannly we CAN also simulate 
different correlation patterns for the clusters.
The factors for abudnance/correlation scale all proteins in the same way (it is enough for simple simulations to see normalization effects on correlation/cluster patterns).
For more elaborate simulations (gradients, proteins of different fold changes) use scRNA simulations
The number of elements in all <factors> parameters must be equal to the number of clusters (here three). The elements in correlationsSets are itself LISTS. This is since
for ONE CLUSTER we could theoretically have several correlationSets that are effected, or we could have only one set of correlted propteins effected by this cluster
#INIRANGE (what cna we simulate)):
#1.) increasing proteinAbudnance Differences between clsiuters (all proteins scaled with same factor...)
#2.) increasing correlation differences (all corr scaled by same factor/ we can have sets of correlated proteins that ARE and that ARE-NOT cluster specifically scaled)
#3.) increasing cell-abundance differences: e.g. a cluster of rare cells with them making up 10% - 1% of cells, …
(number of clusters is not implemented, it would mean that we also have to increase the number of factors, correlations, ...)
#abundanceFactors=\[1.2, 1.5, 2\]
To use the correlations factors, we also have to set the corrleationSets (below), and define the distributions from which a certain number of proteins is sampled to make up protein sets with certain correlations (drawn from a norm dist given by correlation parameters above)
#correlationFactors=\[1.2,1.3,1.5\]
Correaltion sets is proably always used with the same sets everywhere. But might be interesting to simulate e.g. one correlation the is same in all clusters, and another correlation set which is
cluster specific. The list per cluster can have MAXIMALLY as many entries as there are correlationDists. (So every entry is a list of numbers, every number standing for q set of correlations
, where a set is always ONE distribution in proteinCorrelationDist, with zero beeing the baseline and 1 beeing the first dist declared in proteinCorrelationDist)
#correlationSets=\[\[0,1\],\[1\],\[1\]\]
percentages of the cells in each clsuter, they must sum to 100
If this one is a range, all other percentasges are scaled accordingly (e.g. we simulate a rare population with 10% and 5%,  then the other populations are (70,20),(72.5, 22.5))
BE AWARE: this list is of length additionalClusters + 1
The FIRST number is ALWAYS the baseline clsuter number: so 70% of cells will come from a distirbution that was not cahgned by any factors/ with correlations that were not cahnged and come from the unchanged correlation distirbutions
#cellPercentages=\[70,20,10\]
the proteins picked for numberProteins and the proteins picked for correlations (in third parameter of the proteinCorrelationDist) DO NOT AHVE TO BE THE SAME
all those proteins are picked at random, and could, or could not be the same!!! - therefore repeat simulations may times -
Additionally the proteins DO NOT have to be the same between clusters, again, for every clsuter we choose RANDOM proteins to be scaled by the abundanceFactor
#numberProteins=\[2,5,10\]
additionalClusters is the number of additional clsuters next to 'the main' one. So if additionalClusters is equal to one, we have 2 clusters, where the second one is created by
using the abundance/ correlation factors to make this cluster look different to the first one
#additionalClusters=3

## TODO:
- size of neg. binom shouold maybe also vary
- implement bimodal protein counts
- treatmenteffects: we could have in baseline state two clusters, and both r effected by a treatent (maybe even in different ways)



** IMPLEMENTED NORM METHODS**
- TMM
- LEAVE_ONE_OUT_TMM
- CLR_COMPOSITIONS
- CLR_SEURAT
- SUBSAMPING
- LEAVE_ONE_OUT_SUBSAMPLING
- LIBSIZE
- SCTRANSFORM
- SCVI ( scVI from scVI tools, slightly modified code from https://github.com/jmbreda/Sanity/tree/master/reproducibility Sanity benchmark)
- SANITY (bayesian model Sanity: cpp implementation and needs some installations)

implement Sanity:
compile cpp lib// then run python in bhtsne

## ADDITIONAL REQUIREMENTS

**Sanity Normalisation**
for now sanity is not available in e.g. Seurat, and we need to download it from git, compile and adapt python wrappers that r provided in the github-repo

1.) Go to (./src/methods/SanityNormalization)
    git clone https://github.com/jmbreda/Sanity.git
    cd Sanity/src
    make
    (for mac, adapt makefile and add -Xpreprocessing & -lomp for include/linking)
    



## SIMULATIONS

**HOW TO**:
1.) source ./scRNAseq/activate
2.) to run a whole simulation run the script: ./src/simulation/ABseqSimulationBenchmark/simBenchmark.py

PARAMETERS:
    parser.add_argument('--t',   help='threads',
     type=int,   default=-1)
    parser.add_argument('dir',   metavar='DIR',   type=str)
    parser.add_argument('-s',   action='store_false')
    parser.add_argument('-k',   action='store_true')
    parser.add_argument('--d',   help='duplicates',
    
    ### RUN A BENCHMARK WITH INIRANGE INI FILES
    python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./test/data/SimulationFolder/test_simulation.ini --d 2 -k
    
    ### Simulate only groundtruth and simualtions(noise, libsize on top)
        python3 ./src/simulation/main.py ../datasets/Simulations/test.ini --d 1 -k
    
    ### RUN A SINGLE SIMULATION WITH ONE SIUNLATION FOLDER: -s runs on a folder, the folder cna also contain only one file
    (for now running on a file only is not working)
        python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./test/data/SimulationFolder/test_simulation.ini --d 2 -k -s

    ### Benchmark a set of normalized files
    (you need a directory with normalised files)
    python3 src/benchmark/main.py --groundtruth --iniFile <origional inifile to parse true correlations> --stdout <tdout> --t -1 <folder with norm data tables>

WORKFLOW:

When running a Benchmark: we MUST have a INIRANGE in ini file

1.) AbseqSimulationBenchmark/main:
    in MAIN method we create the TmpDir and inis for every simulation, then call run for every ini file
- creates within dir of simulation.ini a folder TmpIniDIR with seperate inis for every simulation (simulatio.ini contains rules how to create many ini for MAX ONE variable with INIRANGE that states whih variable varies for simulation)
    - creats in BENCHMARKED_DATASETS a folder with current data and the name of the master ini, in this case simulation for all the results of the benchmark

from main we cann RUNSIMULATION:
   creates Benchmark object with run method: runs simualtions, normalisation, and copying data around
   
   RUN method starts data simulations and subsewuent normalisation
    - first step: call simulations in Benchmark.run() 
2.) SIMULATION
  - it runs the AbSimulation/main:
  -   creates Simulation/ Groundtruth file in the folder bin/SIMULATION

3.) NORMALIZATION


4.) Benchmarking

from within RUN method becnhmark is called, which calculated all kinds of metrics on the normalized data

FINAL DATA:
we have one folder for every value in inirange and the dublicate (pattern: <INIVALUE>_<DUBLICATEID>)
finally a few result files from every individual simulation are copied to the final RESULT file, where a column SimulationID is added, which is the name of the fodlers,
so a combination of INIVALUE_DUBLICATE

**SAMPLING NEG BINOM and CORRELATIONS**
- Chelovski method to generated correlTED VALUES
- Order values from gaussioan dist + sample values from neg. binom and assign same rank as gaussian to create neg.Binom from correlated samples
