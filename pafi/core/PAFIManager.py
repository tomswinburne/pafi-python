import itertools
import numpy as np
from mpi4py import MPI
from typing import TypeVar, Generic, Any
from .Parser import Parser
from .ResultsHolder import ResultsHolder
from .BaseManager import BaseManager
from .PAFIWorker import PAFIWorker
from .BaseGatherer import BaseGatherer

class PAFIManager(BaseManager):
    def __init__(self, world: MPI.Intracomm, xml_path: str ) -> None:
        super().__init__(world, xml_path, PAFIWorker, BaseGatherer)
    
    def run(self)->None:
        if self.rank==0:
            print(f"""
            Initialized {self.nWorkers} workers with {self.CoresPerWorker} cores
            <> == time averages,  av/err over ensemble
            """)
        print("WR",self.worker_rank)
        for axes_coord in itertools.product(*self.params.axes.values()):
            dict_axes = dict(zip(self.params.axes.keys(), axes_coord))
            #print(dict_axes)
            results = ResultsHolder()
            results.set_dict(dict_axes)

            if results("Temperature")<0.1:
                # Useful helper for including zero temperature cheaply...
                results.set("SampleSteps",50)
                results.set("ThermSteps",10)
                results.set("ThermWindow",10)
            

            final_results = self.Worker.sample(results)
            self.world.Barrier()
            res = final_results.get_dict(["ReactionCoordinate","aveF"])
            print(self.rank,res)
        
    