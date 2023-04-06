import numpy as np
import os
from mpi4py import MPI
from typing import TypeVar, Generic, Any, List
from ..parsers.Parser import Parser

class BaseGatherer:
    def __init__(self,params:Parser,nWorkers:int,suffix:int,rank:int)->None:
        pass