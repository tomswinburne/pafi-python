import numpy as np
import os
from mpi4py import MPI
from typing import TypeVar, Generic, Any, List
from .Parser import Parser
from scipy.interpolate import CubicSpline

class BaseWorker:
    """
        Data Gatherer for PAFI
    """
    def __init__(self,comm : MPI.Intracomm,params:Parser,tag:int) -> None:
        self.tag = tag
        self.comm = comm
        self.local_rank = comm.Get_rank()
        self.params = params
        self.error_count = 0
        self.scale = np.ones(3)
        self.out_width=16
        self.natoms=0
        self.nlocal=0
        self.offset=0
        self.x=None
        self.name="BaseWorker"
        self.has_errors = False
    
    def load_config(self,file_path:str) -> np.ndarray:
        """
            load in file and return configuration
            Will be LAMMPS....
        """
        assert os.path.exists(file_path)
        return np.loadtxt(file_path)

    def pbc(self,X:np.ndarray)->np.ndarray:
        """
            Minimum image of X. Requires cell...
        """
        return X
    def pbc_dist(self,X:np.ndarray,axis:None|int=None)->float|np.ndarray:
        return np.linalg.norm(self.pbc(X),axis=axis)
    
    def make_path(self):
        """
            TODO: parallel i/o and splining
            Only really a problem with memory limitations, say 2GB / core.
            This implies 300M coordinates at double precision == 1M atoms, 100 planes.
            Typical large-scale use - 150k atoms, 20 planes.
            So memory-heavy but nothing problematic so far, leaving for future
        """
        
        # load configurations
        pc = self.params.PathwayConfigurations
        all_X = [self.load_config(p) for p in pc]

        # determine distance TODO: symmetric??
        if self.params("RealMEPDist"):
            self.r_dist = np.array([self.pbc_dist(X-all_X[0])] for X in all_X)
            self.r_dist /= r_dist[-1]
        else:
            self.r_dist = np.linspace(0.,1.,all_X.shape[0])

        # splining - thank you scipy for 'axis' !!
        all_X = np.array([X.flatten() for X in all_X]) # shape=(nknots,natoms)
        bc = self.params("CubicSplineBoundaryConditions")
        assert bc in ['clamped','not-a-knot','natural']
        self.Spline_X = CubicSpline(self.r_dist,all_X,axis=0,bc_type=bc)
        del all_X # save a bit of memory

    def pathway(self,r:float,nu:int=0,
                scale:float|np.ndarray[3]=1.0)->np.ndarray:
        """
            evaluate the pathway
        """
        if isinstance(scale,float):
            scale = np.eye(3)*float
        else:
            scale = np.diag(scale)
        r = max(r,self.r_dist.min())
        r = min(r,self.r_dist.max())
        return self.Spline_X(r,nu=nu).reshape((-1,3))@scale

    def close(self)->None:
        pass
            


        

