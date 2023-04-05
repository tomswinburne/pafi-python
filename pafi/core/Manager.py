import numpy as np
from mpi4py import MPI
from .Parser import Parameters

class Manager:
    """
        The global PAFI Manager. 
    """
    def __init__(self,world : MPI.Intracomm, xml_path:str) -> None:
        self.rank = world.Get_rank()
        self.nProcs = world.Get_size()
        # Read in configuration file
        self.parser = Parameters(xml_path=xml_path)
        self.CoresPerWorker = int(self.parser("CoresPerWorker"))
        if self.nProcs%self.CoresPerWorker!=0:
            if self.rank==0:
                print(f"""
                    CoresPerWorker={self.CoresPerWorker} must factorize nProcs={self.nProcs}!!
                """)
                exit(-1)
        
        # Establish Workers
        self.nWorkers = self.nProcs // self.CoresPerWorker
        self.instance = self.rank // self.CoresPerWorker
        self.local_rank = self.rank % self.CoresPerWorker
        self.min_valid = int(self.parser("RedoThresh")*self.nWorkers*self.parser("nRepeats"))
        
        # Create worker communicator 
        self.instance_comm = world.Split(self.instance,0)

        # Create global communicator between local_rank == 0 processes
        master_ranks = [i*self.CoresPerWorker for i in range(self.nWorkers)]
        self.ensemble_comm = world.Create(world.group.Incl(master_ranks))

        # seed each worker
        self.parser.seed(self.instance)

        # set up Gatherer to accumulate data 






        
        pass