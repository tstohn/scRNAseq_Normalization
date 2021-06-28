from scipy.stats import nbinom

#data from paper
    # 25 cycles of PCR

    #check average number of proteins in single cell with mean and std of distribution
    #check number of UMIs detected
    #size 3.83 && mu 4.68
    #most Abcounts are between 100 and 500, with a few beeing greater, and IGTB1 beeing around 23000

    #remove a fraction bcs of PCR amplification/ sequencing error

#make ABs anneal to proteins in cell
#protein amount is neg. binomially distributed

#loose a fraction by washing

#Amplify the UMIs by PCR randomly

#sequence a random number of those UMI copies

#introduce a biological variance matrix: e.g. neg. binomial distributed
#introduce a technical variance matrix: AB amount differences, annealing differences
class AntiBodyRange():
    absNumber = None
    mean = None
    std = None

    abCountMatrix = None

    def __init__(self, number, mean, std):
        for i in range (number):
            abCountVector = nbinom.rvs(mean ,std, size = number)



class Parameters():

    #incubation time is considered constant for now
    #antibody count is not considered for now: in sc IDseq constant AB concentration was used for all targets (0.1myg/ml over night)
    #biological variation matrix
    #technical variation matrix
    """ Raw Data Parameters """
    proteinCountMatrix = None
    """ Experimental Parameters """
    abBindingEfficiency = None
    pcrCycles = None
    seqEfficiency = None

    """ PARAMETERS
        rangeVector: a vector of quadruples(range, number, mean, std) for Abs of different distributions
        abDuplicates: to simulate the usage of several tags for one AB
    """

    def __simulateProteinCountMatrix(self, rangeVector, abDuplicates = 1):
        for range in rangeVector:
            abRange = AntiBodyRange(range)

    def __init__(self, rangeVector):
        self.__simulateProteinCountMatrix(range)

class SingleCellSimulation():

    parameters = None
    #sample * AB matrix of simulated data
    simulatedData = None

    def __init__(self, parameters):
        self.parameters = parameters

    def simulateData(self):
        #simulate AB binding efficiency

        #simulate PCR amplification and sequencing
        return(self.simulatedData)
