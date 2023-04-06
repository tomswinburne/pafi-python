import numpy as np
import os
from typing import TypeVar, Generic, Any, List
from .Parser import Parser
from mpi4py import MPI
from lammps import lammps, MPIAbortException

from .BaseWorker import BaseWorker

class LAMMPSWorker(BaseWorker):
    def __init__(self, 
                 comm: MPI.Intracomm, params: Parser, tag: int) -> None:
        super().__init__(comm, params, tag)
        """
            Start up LAMMPS worker
        """
        self.name = "LAMMPSWorker"
        self.last_error_message = ""
        
        self.load_lammps()        
        if self.has_errors:
            return
        
        self.reset()
        start_config = self.params.PathwayConfigurations[0]
        self.run_script("Input",{'FirstPathImage',start_config})
        if self.has_errors:
            return
        
        self.has_cell_data = False
        self.lammps_setup()
        
        self.fill_lammps_vectors()
    
    def load_lammps(self)->None:
        """
            Load LAMMPS instance
        """
        if self.params("LogLammps"):
            logfile = 'log.lammps.%d' % self.tag 
        else:
            logfile = 'none'
        try:
            cmdargs = ['-screen','none','-log',logfile]
            self.L = lammps(comm=self.comm,cmdargs=cmdargs)
            self.check_lammps_compatibility()
        except MPIAbortException as ae:
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
    
    def reset(self)->None:
        """
            reset the box
        """
        self.run_commands("clear")
        self.scale = np.ones(3)
        self.made_fix=False
        self.made_compute=False
    
    def run_script(self,key:str,args:dict)->None:
        """
            Run a script defined in the XML
            Important to replace any %wildcards% if they are there!
            TODO - catch these?
        """
        if key in self.params.scripts:
            script = self.params.parse_script(key,args)
            self.run_commands(script)

    def run_commands(self,cmds : str | List[str]) -> bool:
        """
            Run LAMMPS commands line by line, checking for errors
        """
        cmd_list = cmds.splitlines() if isinstance(cmds,str) else cmds
        for cmd in cmd_list:
            try:
                self.L.command(cmd)
            except MPIAbortException as ae:
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
        except MPIAbortException as ae:
            if self.local_rank==0:
                print("Error in gather:",ae)
            self.last_error_message = ae
        return res
    
    def scatter(self,name:str,data:np.ndarray):
        """
            Scatter data to LAMMPS
            Assume ordered with ID
        """
        if np.issubdtype(data.dtype,int):
            type = 0
        elif np.issubdtype(data.dtype,float):
            type = 1
        assert len(data.shape)==2
        count = data.shape[1] if len(data.shape)>1 else 1
        try:
            self.L.scatter(name,type,count,data.flatten())
        except MPIAbortException as ae:
            if self.local_rank==0:
                print("Error in scatter:",ae)
            self.last_error_message = ae
    
    def lammps_setup(self)->None:
        """
            Extract natoms, cell etc.
        """
        self.natoms = self.L.get_natoms()
        # boxlo, boxhi, xy, yz, xz, periodicity, box_change
        # TODO periodicity flags?
        cell_data = self.L.extract_box()
        self.Periodicity = np.r_[cell_data[5]].astype(bool)
        self.Cell = np.diag(np.r_[[cell_data[1]]]-np.r_[[cell_data[0]]])
        self.Cell[0][1] = cell_data[2]
        self.Cell[1][2] = cell_data[3]
        self.Cell[0][2] = cell_data[4]
        self.invCell = np.linalg.inv(self.Cell)
        self.has_cell_data = True
    
    def pbc(self,X:np.ndarray,central:bool=True)->np.ndarray:
        """
            Minimum image convention, using cell data
            central : bool
                map scaled coordinates to [-.5,.5] if True, else [0,1]
        """
        if not self.has_cell_data:
            super().pbc(X)
        else:
            sX = X.reshape((-1,3))@self.invCell
            sX -= np.floor(sX+0.5*central)@np.diag(self.Periodicity)
            return (sX@self.Cell).reshape((X.shape))

    def load_config(self, file_path: str) -> np.ndarray:
        self.made_fix = False
        self.make_compute = False
        self.run_commands(f"""
            delete_atoms group all
            read_data {file_path} add merge
        """)
        return self.gather("x").reshape((-1,3))

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

    def populate(self,r:float,T:float)->float:
        """
            Establish worker on a given plane with the reference pathway
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
        norm_t = np.linalg.norm(path_t)
        path_t /= norm_t
        for i,c in enumerate(["d_nx","d_ny","d_nz"]):
            self.scatter(c,path_t[:,i])
        del path_t

        # fill dtangent: d_dn[x,y,z]
        path_t = self.pathway(r,nu=2,scale=self.scale) / norm_t**2
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

    def extract_compute(self,id:str,vector:bool=True)->Any:
        """
            extract compute from LAMMPS
            id : str
                compute id
            vector: bool
                is the return value a vector

        """
        from lammps import LMP_STYLE_GLOBAL,LMP_TYPE_VECTOR,LMP_TYPE_SCALAR
        style = LMP_STYLE_GLOBAL
        type = LMP_TYPE_VECTOR if vector else LMP_TYPE_SCALAR
        """
            Wrapper of LAMMPS extract_compute. 
            Using `import lammps` to prevent clashes with mulitple installations
        """
        assert hasattr(self.L,"numpy")
        try:
            res = self.L.numpy.extract_compute(id,style,type) # type: ignore
            return np.array(res)
        except Exception as e:
            if self.rank==0 and self.verbose:
                print("FAIL EXTRACT COMPUTE",e)
            self.close()
            return None
    
    def getEnergy(self)->float:
        """
            Extract the potential energy
        """
        return self.extract_compute("thermo_pe",vector=False)


    def close(self) -> None:
        """
            Close down. TODO Memory management??
        """
        super().close()
        self.L.close()
        