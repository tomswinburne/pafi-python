from ..parsers.Parser import Parser

class BaseGatherer:
    def __init__(self,params:Parser,nWorkers:int,rank:int)->None:
        """Basic gatherer of PAFI simulation data

        Parameters
        ----------
        params : Parser
            Custom or predefined PAFI Parser object
        nWorkers : int
            total number of PAFI workers
        rank : int
            global MPI rank
        """
        self.params = params
        self.nWorkers = nWorkers
        self.rank = rank