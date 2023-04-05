import numpy as np
from mpi4py import MPI

class GenericWorker:
    """
        Generic PAFI worker, with all non-LAMMPS functionality
    """

    def __init__(self,instance_comm,parser,)