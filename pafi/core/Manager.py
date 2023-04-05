import numpy as np
from mpi4py import MPI


class PAFIManager:
    """
        The global PAFI Manager. 
    """
    def __init__(self,comm : MPI.Intracomm,config_file:str) -> None:
        rank = comm.Get_rank()
        nProcs = comm.Get_size()
        pass