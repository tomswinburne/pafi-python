from mpi4py import MPI
import os
from typing import List
from ..parsers.Parser import Parser
from .ResultsHolder import ResultsHolder

class BaseGatherer:
    def __init__(self,params:Parser,
                 nWorkers:int,rank:int,
                 ensemble_comm:MPI.Intracomm,
                 roots:List[int])->None:
        """Basic gatherer of PAFI simulation data

        Parameters
        ----------
        params : Parser
            Custom or predefined PAFI Parser object
        nWorkers : int
            total number of PAFI workers
        rank : int
            global MPI rank
        ensemble_comm : MPI.Intracomm
            MPI communicator to gather ensemble data
        roots : List[int]
            list of root ranks for each worker. Only collate data here
            This should be depreceated asap
        """
        self.params = params
        self.nWorkers = nWorkers
        self.rank = rank
        self.comm = ensemble_comm
        self.roots = roots
        self.epoch_data = None # for each cycle
        self.all_data = None # for total simulation
    
    def gather(self,data:dict|ResultsHolder)->None:
        """Gather results from a simulation epoch,
        local to each worker. Here, very simple,
        is overwritten by each call of gather()

        Parameters
        ----------
        data : dict or ResultsHolder
            Simulation data, extracted as dictionary from ResultsHolder
        """
        if isinstance(data,ResultsHolder):
            self.epoch_data = data.data.copy()
        else:
            self.epoch_data = data.copy()
       
    
    def collate(self):
        """Collate all data on root node
        It is assumed that all data is in dictionary form, with identical keys
        Only basic checks for multiple calls
        """
        if not self.epoch_data is None:
            all_epoch_data = self.comm.gather(self.epoch_data)
            self.epoch_data = None
        else:
            all_epoch_data = None

        if self.rank == 0 and not all_epoch_data is None:
            if self.all_data is None:
                self.all_data = {k:[] for k in all_epoch_data[0].keys()}
            for k in all_epoch_data[0].keys():
                self.all_data[k] += list(d[k] for d in all_epoch_data)
    
    def write_pandas(self,path:os.PathLike[str])->None:
        """Write data as pandas dataframe

        Parameters
        ----------
        path : os.PathLike[str]
            path to file
        """

        if self.rank==0:
            if self.all_data is None:
                print("No data to write! Exiting!")
            else:
                import pandas as pd
                if not os.path.isdir(os.path.dirname(path)):
                    raise IOError("Unknown directory for writing csv!")
                else:
                    df = pd.DataFrame(self.all_data)
                    df.to_csv(path)
            
                
