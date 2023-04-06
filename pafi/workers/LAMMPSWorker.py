import numpy as np
import os
from typing import Any, List
from ..parsers.Parser import Parser
from mpi4py import MPI
from lammps import lammps,LMP_STYLE_GLOBAL,LMP_TYPE_VECTOR,LMP_TYPE_SCALAR
from .BaseWorker import BaseWorker
from ..results.ResultsHolder import ResultsHolder

class LAMMPSWorker(BaseWorker):
    """LAMMPS worker for PAFI, inheriting BaseWorker

        Initialization routine:
        - runs "Input" script, loading first configuration
        - extracts cell data
        - defines various functions wrapping computes etc.
        - defines initialize_hyperplane function to set plane

        Parameters
        ----------
        comm : MPI.Intracomm
            MPI communicator
        params : Parser
            Predefined or custom PAFI Parser object
        worker_instance : int
            unique worker rank
        """
    def __init__(self, comm: MPI.Intracomm, 
                 params: Parser, worker_instance: int) -> None:
        super().__init__(comm, params, worker_instance)
        self.name = "LAMMPSWorker"
        self.last_error_message = ""
        self.start_lammps()
        if self.has_errors:
            print("ERROR STARTING LAMMPS!",self.last_error_message)
            return
        start_config = self.params.PathwayConfigurations[0]
        # TODO abstract
        self.run_script("Input")
        if self.has_errors:
            print("ERROR RUNNING INPUT SCRIPT!",self.last_error_message)
            return
        
        self.has_cell_data = False
        self.get_cell_data()
        # See initialize_hyperplane
        self.scale = np.ones(3)
        self.made_fix=False
        self.made_compute=False
        self.make_path()
        
    
    def start_lammps(self)->None:
        """Initialize LAMMPS instance

            Optionally 

        """
        if self.params("LogLammps"):
            logfile = 'log.lammps.%d' % self.worker_instance 
        else:
            logfile = 'none'
        try:
            cmdargs = ['-screen','none','-log',logfile]
            self.L = lammps(comm=self.comm,cmdargs=cmdargs)
            self.check_lammps_compatibility()
        except Exception as ae:
            print("Couldn't load LAMMPS!")
            self.has_errors = True
    
    def check_lammps_compatibility(self)->None:
        """
            Ensure LAMMPS is new enough and has fix_pafi
        """
        lammps_release_int = self.L.version()
        if lammps_release_int<20201101: 
            if self.local_rank == 0:
                print("Require LAMMPS version > 01Nov2020!")
            return
        pafi_package = "EXTRA-FIX" if lammps_release_int>=20210728 else "USER-MISC"
        self.has_pafi = self.L.has_package(pafi_package)
        if not self.has_pafi and self.local_rank==0:
            print("Cannot find PAFI package in LAMMPS!")
            self.has_errors = True
        
    def run_script(self,key:str,arguments:None|dict=None)->None:
        """
            Run a script defined in the XML
            Important to replace any %wildcards% if they are there!
            TODO - catch these?
        """
        if key in self.params.scripts:
            script = self.params.parse_script(key,arguments=None)
            self.run_commands(script)

    def run_commands(self,cmds : str | List[str]) -> bool:
        """
            Run LAMMPS commands line by line, checking for errors
        """
        cmd_list = cmds.splitlines() if isinstance(cmds,str) else cmds
        for cmd in cmd_list:
            try:
                self.L.command(cmd)
            except Exception as ae:
                if self.local_rank==0:
                    print("LAMMPS ERROR:",cmd,ae)
                self.last_error_message = ae
    
    def gather(self,name:str,type:None|int=None,count:None|int=None):
        """
            Wrapper of LAMMPS gather()
            will be ordered by ID
        """
        if name in ['x','f','v'] or 'f_' in name:
            if type is None:
                type = 1
            if count is None:
                count = 3
        elif name in ['id','type','image']:
            if type is None:
                type = 0
            if count is None:
                count = 1
        if type is None or count is None:
            raise ValueError("Error in gather: type or count is None")
        
        try:
            res = self.L.gather(name,type,count)
        except Exception as ae:
            if self.local_rank==0:
                print("Error in gather:",ae)
            self.last_error_message = ae
        return np.ctypeslib.as_array(res).reshape((-1,count))
    
    def scatter(self,name:str,data:np.ndarray):
        """
            Scatter data to LAMMPS
            Assume ordered with ID
        """
        if np.issubdtype(data.dtype,int):
            type = 0
        elif np.issubdtype(data.dtype,float):
            type = 1
        count = data.shape[1] if len(data.shape)>1 else 1
        try:
            self.L.scatter(name,type,count,
                           np.ctypeslib.as_ctypes(data.flatten()))
        except Exception as ae:
            if self.local_rank==0:
                print("Error in scatter:",ae)
            self.last_error_message = ae
    
    def get_natoms(self)->int:
         """
            Get the atom count
         """
         return self.L.get_natoms()
    
    def get_cell_data(self)->None:
        """
            Extract supercell
        """
        boxlo,boxhi,xy,yz,xz,pbc,box_change = self.L.extract_box()
        self.Periodicity = np.array([bool(pbc[i]) for i in range(3)],bool)
        self.Cell = np.zeros((3,3))
        for cell_j in range(3):
            self.Cell[cell_j][cell_j] = boxhi[cell_j]-boxlo[cell_j]
        self.Cell[0][1] = xy
        self.Cell[0][2] = xz
        self.Cell[1][2] = yz
        self.invCell = np.linalg.inv(self.Cell)
        self.has_cell_data = True

    def load_config(self, file_path: str) -> np.ndarray:
        self.made_fix = False
        self.make_compute = False
        self.run_commands(f"""
            delete_atoms group all
            read_data {file_path} add merge
        """)
        return self.gather("x",type=1,count=3)

    def thermal_expansion_supercell(self,T:float=0) -> None:
        """
            rescale the LAMMPS supercell NOT COORDINATTES
            according to the provided thermal expansion data
        """
        newscale = self.params.expansion(T)
        rs = newscale / self.scale
        self.run_commands(f"""
            change_box all x scale {rs[0]} y scale {rs[1]} z scale {rs[2]}
            run 0""")
        self.scale = newscale.copy()

    def extract_compute(self,id:str,vector:bool=True)->Any:
        """
            extract compute from LAMMPS
            id : str
                compute id
            vector: bool
                is the return value a vector

        """
        style = LMP_STYLE_GLOBAL
        type = LMP_TYPE_VECTOR if vector else LMP_TYPE_SCALAR
        assert hasattr(self.L,"numpy")
        try:
            res = self.L.numpy.extract_compute(id,style,type) 
            return np.array(res)
        except Exception as e:
            if self.local_rank==0:
                print("FAIL EXTRACT COMPUTE",e)
            self.close()
            return None
    
    def extract_fix(self,id:str,size:int=1)->Any:
        """
            extract fix from LAMMPS
            id : str
                fix id
            index: int
                return index if vector
        """
        style = LMP_STYLE_GLOBAL
        type = LMP_TYPE_VECTOR if size>1 else LMP_TYPE_SCALAR
        assert hasattr(self.L,"numpy")
        try:

            res = lambda i: self.L.numpy.extract_fix(id,style,type,ncol=i)
            if size>1:
                return np.array([res(i) for i in range(size)])
            else:
                return res(0)
        except Exception as e:
            if self.local_rank==0:
                print("FAIL EXTRACT FIX",id,e)
            self.close()
            return None
    
    def get_energy(self)->float:
        """
            Extract the potential energy
        """
        return self.extract_compute("thermo_pe",vector=False)

    def initialize_hyperplane(self,r:float,T:float)->None:
        """
            Establish worker on a given plane with the 
            reference pathway. Returns the tangent magnitude
        """
        self.thermal_expansion_supercell(T) # updates self.scale also

        # check for __pafipath fix to store the path data
        # (path : d_u[x,y,z],tangent: d_n[x,y,z],dtangent: d_dn[x,y,z])
        if not self.made_fix:
            self.run_commands(f"""
            fix __pafipath all property/atom d_ux d_uy d_uz d_nx d_ny d_nz d_dnx d_dny d_dnz
            run 0""")
            self.made_fix=True
        
        # fill positions: x and d_u[x,y,z]
        path_x = self.pathway(r,nu=0,scale=self.scale)
        self.scatter("x",path_x)
        for i,c in enumerate(["d_ux","d_uy","d_uz"]):
            self.scatter(c,path_x[:,i])
        del path_x

        # fill tangent: d_n[x,y,z]
        path_t = self.pathway(r,nu=1,scale=self.scale)
        path_t -= path_t.mean(0)
        self.norm_t = np.linalg.norm(path_t)
        path_t /= self.norm_t
        for i,c in enumerate(["d_nx","d_ny","d_nz"]):
            self.scatter(c,path_t[:,i])
        del path_t

        # fill dtangent: d_dn[x,y,z]
        path_t = self.pathway(r,nu=2,scale=self.scale) 
        
        path_t /= self.norm_t**2
        for i,c in enumerate(["d_dnx","d_dny","d_dnz"]):
            self.scatter(c,path_t[:,i])
        del path_t

        self.run_commands("run 0")

        # check for __pafipath compute to make this data accessible
        if not self.made_compute:
            self.run_commands(f"""
            compute __pafipath all property/atom d_ux d_uy d_uz d_nx d_ny d_nz d_dnx d_dny d_dnz
            run 0""")
            self.made_compute=True
    
    def close(self) -> None:
        """
            Close down. TODO Memory management??
        """
        super().close()
        self.L.close()

    