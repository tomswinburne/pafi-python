import numpy as np
import pandas as pd
from mpi4py import MPI
from typing import Any, List
from ..parsers.Parser import Parser
from .BaseGatherer import BaseGatherer
class Gatherer(BaseGatherer):
    """Gatherer of PAFI simulation data

        Parameters
        ----------
        params : Parser
            Custom or predefined PAFI Parser object
        nWorkers : int
            total number of PAFI workers
        rank : int
            global MPI rank
        """
    def __init__(self, params: Parser, nWorkers: int, rank: int) -> None:
        super().__init__(params, nWorkers, rank)
