class Cluster:
    def __init__(self, data):
        self.data = data

#represent cells by a function of cell state (cell stae is a continuous space)
#model this space as a manifold, create a manufold for each treatment and make a mapping of 
#a manifold distriution function describing all possible cell states to this function for a different treatment

#put those mappings into a auto-encoder to elarn the hidden continuous space of cell state points for different treatments
class CellClustering:
    def __init__(self, data):
        self.data = data

    def meld_algorithm(self):
        
