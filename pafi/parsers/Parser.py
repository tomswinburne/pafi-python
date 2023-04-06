import numpy as np
import os,glob
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
from typing import Union,Any
ScriptArg = Union[int,float,str]

from .BaseParser import BaseParser

class Parser(BaseParser):
    """Default reader of PAFI XML configuration file

        Parameters
        ----------
        xml_path : str
            path to XML file
        
        Methods
        ----------
        __call__
        find_suffix_and_write
        

        Raises
        ------
        IOError
            If path is not found
    """
    def __init__(self, xml_path: str) -> None:
        super().__init__(xml_path)
        
        # initial seed, but must be different across workers...
        self.seeded = False
        
        # repeats, valid bounds (see set_min_valid())
        self.maxRepeats = self.parameters["nRepeats"]
        self.maxRepeats += self.parameters["maxExtraRepeats"]         
        self.set_min_valid(1) # temporary
    
    def set_min_valid(self,nWorkers):
        self.minValidResults = self.parameters["ReSampleThresh"]
        self.minValidResults *= self.parameters["nRepeats"]
        self.minValidResults *= nWorkers
        self.minValidResults = int(self.minValidResults)
    
    def min_valid(self):
        return self.minValidResults
    
    def seed(self,worker_instance:int)->None:
        """Generate random number seed

        Parameters
        ----------
        worker_instance : int
            unique to each worker, to ensure independent seed
        """
        if not self.seeded:
            self.randseed = self.parameters["GlobalSeed"] * (worker_instance+1)
            self.rng = np.random.default_rng(self.randseed)
            self.rng_int = self.rng.integers(low=100, high=10000)
            self.seeded=True
    
    def randint(self)->int:
        """
            Generate random integer.
            Gives exactly the same result each time unless reseed=True
        """
        if not self.seeded:
            print("NOT SEEDED!!")
            exit(-1)
        
        if self.parameters["FreshSeed"]:
            self.rng_int = self.rng.integers(low=100, high=10000)
        return str(self.rng_int)
    
    def expansion(self,T:float)->np.ndarray:
        """Return the relative x,y,z thermal expansion,
            using the data provided in the XML file

        Parameters
        ----------
        T : float
            Temperature in K

        Returns
        -------
        np.ndarray, shape (3,)
            relative x,y,z thermal expansion
        """
        scale = np.ones(3) 
        scale += self.parameters["LinearThermalExpansion"]*T
        scale += self.parameters["QuadraticThermalExpansion"]*T*T
        return scale
        

            


