import itertools
import numpy as np
from mpi4py import MPI
from typing import TypeVar, Generic, Any
from .Parser import Parser
from .ResultsHolder import ResultsHolder
from .BaseWorker import BaseWorker
from .BaseGatherer import BaseGatherer

class BaseManager:
    """
        The global PAFI Manager. 

        run :
    """
    def __init__(self,world : MPI.Intracomm, xml_path:str,
                 Worker=BaseWorker,Gatherer=BaseGatherer) -> None:
        self.world = world
        self.rank = world.Get_rank()
        self.nProcs = world.Get_size()
        # Read in configuration file
        self.params = Parser(xml_path=xml_path)
        self.CoresPerWorker = int(self.params("CoresPerWorker"))
        if self.nProcs%self.CoresPerWorker!=0:
            if self.rank==0:
                print(f"""
                    CoresPerWorker={self.CoresPerWorker} must factorize nProcs={self.nProcs}!!
                """)
                exit(-1)
        """
            Establish Workers
            worker_comm : Worker communicator for e.g. LAMMPS
        """
        self.nWorkers = self.nProcs // self.CoresPerWorker
        
        # Create worker communicator 
        self.worker_rank = self.rank // self.CoresPerWorker
        self.worker_comm = world.Split(self.worker_rank,0)
        
        # set up and seed each worker
        self.params.seed(self.worker_rank)
        self.Worker = Worker(self.worker_comm,
                             self.params,
                             self.worker_rank)
        
        """
            Establish Gatherer
            ensemble_comm: Global communicator for averaging
        """
        self.roots = [i*self.CoresPerWorker for i in range(self.nWorkers)]
        self.ensemble_comm = world.Create(world.group.Incl(self.roots))
        

        self.Gatherer = None
        if self.rank in self.roots:
            worker_errors = self.ensemble_comm.gather(self.Worker.has_errors)
            if self.rank==0 and max(worker_errors)>0:
                raise IOError("Worker Errors!")
            self.Gatherer = Gatherer(self.params,
                                 self.nWorkers,
                                 self.params.suffix,
                                 self.rank)
        if self.rank==0:
            print(self.params.welcome_message())   
    
    def close(self)->None:
        self.Worker.close()
            