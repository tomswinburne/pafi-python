import numpy as np
import os,glob
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
from typing import Union,Any
ScriptArg = Union[int,float,str]

from .BaseParser import BaseParser

class Parser(BaseParser):
    """
        BaseParameters for a PAFI run,
        read from an XML file.
        with a few extra helper functions
    """
    def __init__(self, xml_path: str) -> None:
        super().__init__(xml_path)
        
        # initial seed, but must be different across workers...
        self.rng = np.random.default_rng(self.parameters['GlobalSeed'])
        self.seeded = False
        
        # repeats, valid bounds (see set_min_valid())
        self.maxRepeats = self.parameters["nRepeats"]
        self.maxRepeats += self.parameters["maxExtraRepeats"]         
        self.set_min_valid(1) # temporary
    
    def set_min_valid(self,nWorkers):
        self.minValidResults = self.parameters["RedoThresh"]
        self.minValidResults *= self.parameters["nRepeats"]
        self.minValidResults *= nWorkers
        self.minValidResults = int(self.minValidResults)
    
    def min_valid(self):
        return self.minValidResults
    
    def seed(self,worker_instance):
        """
            Generate random number seed- why is this here...
        """
        self.seed = self.parameters["GlobalSeed"] * (worker_instance+1)
        self.rng = np.random.default_rng(self.seed)
        self.rng_int = self.rng.integers(low=100, high=10000)
        self.seeded=True
    
    def seed_str(self):
        """
            Gives exactly the same seed each time unless reseed=True
        """
        if not self.seeded:
            print("NOT SEEDED!!")
        
        if self.parameters["FreshSeed"]:
            self.rng_int = self.rng.integers(low=100, high=10000)
        return str(self.rng_int)

    def expansion(self,T:float)->np.ndarray:
        """
            Return the relative x,y,z thermal expansion,
            using the data provided in the XML file
        """
        scale = np.ones(3) 
        scale += self.parameters["LinearThermalExpansion"]*T
        scale += self.parameters["QuadraticThermalExpansion"]*T*T
        return scale
        

            


