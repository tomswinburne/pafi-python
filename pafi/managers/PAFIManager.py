import itertools
from typing import List
import numpy as np
import os
from mpi4py import MPI
from ..results.ResultsHolder import ResultsHolder
from .BaseManager import BaseManager
from ..parsers.PAFIParser import PAFIParser
from ..workers.PAFIWorker import PAFIWorker
from ..results.Gatherer import Gatherer

class PAFIManager(BaseManager):
    def __init__(self, world: MPI.Intracomm, 
                 xml_path:None|os.PathLike[str]=None,
                 params:None|PAFIParser=None,
                 Worker:PAFIWorker=PAFIWorker,
                 Gatherer:Gatherer=Gatherer) -> None:
        """Default manager of PAFI, child of BaseManager

        Parameters
        ----------
        world : MPI.Intracomm
            MPI communicator
        xml_path : None or os.PathLike[str]
            path to XML configuration file
        params : None or PAFIParser object
            preloaded PAFIParser object, 
        Worker : PAFIWorker, optional,
            Can be overwritten by child class, by default PAFIWorker
        Gatherer : Gatherer, optional
            Can be overwritten by child class, by default Gatherer
        """
        
        
        assert (not params is None) or (not xml_path is None)
        if params is None:
            params = PAFIParser(xml_path=xml_path)
        super().__init__(world, params, Worker, Gatherer)
        
    
    
    def run(self,print_fields:List[str]|None=None,
            width:int=10,precision:int=5)->None:
        """Basic parallel PAFI sampling

            Performs a nested loop over all <Axes>, in the order
            presented in the XML configuration file.
            Here, parallelization is naive- all workers are given the 
            same parameters.
        Parameters
        ----------
        print_fields : List[str] or None
            Fields to print to screen, default None. 
            If None, will print "Temperature","ReactionCoordinate","aveF"
        width : int
            character count of field printout, default 10
        precision : int
            precision of field printout, default 4
        """
        if print_fields is None:
            print_fields = \
                ["Rank","Temperature","ReactionCoordinate","aveF","vadrF"]
        for f in print_fields:
            width = max(width,len(f))
        
        def line(data,top=False):
            format_string = ""
            format_string = ("{: >%d} "%width)*len(data)
            if isinstance(data,dict):
                _fields = [data[f] for f in print_fields]
            else:
                _fields = data
            
            fields = []
            for f in _fields:
                isstr = isinstance(f,str)
                fields += [f if isstr else np.round(f,precision)]
            return format_string.format(*fields)

        if self.rank==0:
            print(f"""
            Initialized {self.nWorkers} workers with {self.CoresPerWorker} cores
            <> == time averages,  av/err over ensemble
            """)
            print(line(print_fields,top=True))
                
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
            final_results.set("Rank",self.rank)
            if self.rank in self.roots:
                self.Gatherer.gather(final_results)
                self.Gatherer.collate()
            
            

            self.world.Barrier()
            res = final_results.get_dict(print_fields[1:],blanks="n/a")
            res["Rank"] = self.rank
            if self.rank in self.roots:
                all_res_str = self.ensemble_comm.gather(res)
            else:
                all_res_str = None
            self.world.Barrier()
            if self.rank==0:
                for res in all_res_str:
                    print(line(res))
        if self.rank==0:
            pandas_csv = self.params.csv_file
            self.Gatherer.write_pandas(path=pandas_csv)
        
    