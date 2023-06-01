import pandas as pd
import plotext as plt
import os,glob
import numpy as np
from typing import Any,List
from .ResultsHolder import ResultsHolder
from ..parsers.PAFIParser import PAFIParser
from scipy.integrate import cumtrapz

class ResultsProcessor:
    def __init__(self,
                 data_path:os.PathLike[str]|List[os.PathLike[str]],
                 xml_path:None|os.PathLike[str]=None) -> None:
        """Read in PAFI data and plot results

        Parameters
        ----------
        data_path : os.PathLike[str] | List[os.PathLike[str]]
            path, wildcard, or list of paths to PAFI csv files
        xml_path : None | os.PathLike[str]
            path to PAFI XML configuration file, default None
        Methods
        ----------
        integrate()
        """

        # Read in parameters
        self.params = PAFIParser(xml_path)

        # Read in 
        data = []
        if isinstance(data_path,list):
            for dp in data_path:
                if os.path.exists(dp):
                    data += [pd.read_csv(dp)]
        else:
            assert os.path.exists(data_path)
            self.pd_data = pd.read_csv(data_path)

    def integrate(self,
                  argument:str='ReactionCoordinate',
                  target:str='aveF',
                  remesh:int=5)->pd.DataFrame:
        """Cumulative integration of data along an axis. 
        Integration routine makes a spline interpolation
        to increase the number of intergrand evaluations 
        
        Parameters
        ----------
        argument : str, optional
            integration argument, by default 'ReactionCoordinate'
        target : str, optional
            integral, by default 'aveF'
        remesh : int, optional
            the number of integrand evaluations between 
            existing knot points, by default 5
        Returns
        -------
        pd.DataFrame
            return data
        """
        
        x = argument
        y = target

        pd_data = self.pd_data.drop(x,axis=1).drop_duplicates()
        pd_fields = list(pd_data.columns) + [x,y+"_int",y+'_int_std']
        all_data = []
        for p in pd_data.iterrows():
            ee = self.pd_data.copy()
            for key,val in p[1].items():
                ee = ee[ee[key]==val]
            
            data = np.r_[[ee[w].to_numpy() for w in [x,y]]].T
            #data = self.remesh(data[data[:,0].argsort(),:],remesh)
            y_a = np.append(np.zeros(1),cumtrapz(data[:,1],data[:,0]))
            #y_e = np.append(np.zeros(1),cumtrapz(data[:,2],data[:,0]))
            all_data += [[*list(p[1].to_dict().values()),data[:,0],y_a]]#,y_e]]
        return pd.DataFrame(data=all_data,columns=pd_fields)