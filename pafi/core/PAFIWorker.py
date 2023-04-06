import numpy as np
import os
from mpi4py import MPI
from typing import Any, List
from .Parser import Parser
from .LAMMPSWorker import LAMMPSWorker
from .ResultsHolder import ResultsHolder

class PAFIWorker(LAMMPSWorker):
    """
        The hyperplane-constrained sampling run.
        results: ResultsHolder instance
            essentially a dictionary, defining both 
            input and output data- temeperature, reaction coordinate
            and other are defined, and others are filled.
            *Crucially, if a parameter is specified in results, 
            it overrides any value in the parameter file*
            Must have at least `ReactionCoordinate` and `Temperature` defined
        
        The PAFI workflow can be summarised as:
            1) Execute `PreRun` script
            2) Apply `fix_pafi` constraint at defined `ReactionCoordinate`
            3) Execute `PreTherm` script
            4) Thermalization for `ThermSteps` steps at `Temperature`
            5) Execute `constrained_average()` function
                In standard PAFI this runs for `SampleSteps` steps and 
                time averages the output of `fix_pafi`, as shown below.
                See https://docs.lammps.org/fix_pafi.html for details.
            6) Extract the average displacment from path if `PostDump==1`
            7) Minimize in-plane to check system returns to path.
                The check is the max per-atom displacement : `MaxJump`
                If `MaxJump` is larger than `MaxJumpMaxJumpThresh` then 
                the sample is retained but marked as `Valid=False`
            8) Execute `PostRun` script
        
        The `contrained_average()` function is therefore suitable for
        any form of hyperplane-constrained averaging.
        """
    
    def __init__(self, comm: MPI.Intracomm, params: Parser, tag: int) -> None:
        super().__init__(comm, params, tag)
    
    def constrained_average(self,results:ResultsHolder)->ResultsHolder:
        """
            The hyperplane-contrained averaging step. 
            See the PAFIWorker() documentation
            results: ResultsHolder instance
                essentially a dictionary with inputs and outputs
                we read it in, add to it, and return
        """
        # we want any input paramaters in results() to override self.params()
        params = lambda k: results(k) if results.has_key(k) else self.params(k)
        steps = params("SampleSteps")

        # the fix 'pafi' is a fix_pafi instance, constraining the dynamics
        self.run_commands(f"""
            fix avepafi all ave/time 1 {steps} {steps} f_pafi[*]
            run {steps}
        """)

        # helper function to convert avepafi to a dictionary
        pafi_results_dict = self.extract_pafi_data("avepafi")
        
        # add all fields to the results object
        results.set_dict(pafi_results_dict)
        
        # unfix the average
        self.run_commands("unfix avepafi")
        
        return results
    

    """
        Helper functions to perform the sampling

    """
    def extract_pafi_data(self,name:str="avepafi")->dict:
        """
            extract data from fix_pafi and return a dictionary
        """
        res = {}
        fix_data = self.extract_fix(name,size=4)
        res['aveF'] = -fix_data[0] * self.norm_t
        res['varF'] = fix_data[1]**2 * self.norm_t**2 - res['aveF']**2
        res['avePsi'] = fix_data[2]
        res['dXTangent'] = fix_data[3]
        return res
    
    def sample(self,results:ResultsHolder)->ResultsHolder:
        """
        Main sampling run.
        1) Execute `PreRun` script
        2) Apply `fix_pafi` constraint at defined `ReactionCoordinate`
        3) Execute `PreTherm` script
        4) Thermalization for `ThermSteps` steps at `Temperature`
        5) Execute `constrained_average()` function
            In standard PAFI this runs for `SampleSteps` steps and 
            time averages the output of `fix_pafi`, as shown below.
            See https://docs.lammps.org/fix_pafi.html for details.
        6) Extract the average displacment from path if `PostDump==1`
        7) Minimize in-plane to check system returns to path.
            The check is the max per-atom displacement : `MaxJump`
            If `MaxJump` is larger than `MaxJumpMaxJumpThresh` then 
            the sample is retained but marked as `Valid=False`
        8) Execute `PostRun` script
    

        results: ResultsHolder instance
            essentially a dictionary, defining both 
            input and output data- temeperature, reaction coordinate
            and other are defined, and others are filled.
            *Crucially, if a parameter is specified in results, 
            it overrides any value in the parameter file*
            *Must* have `ReactionCoordinate` and `Temperature` defined
        """

        # the esssentials
        assert results.has_key("ReactionCoordinate")
        assert results.has_key("Temperature")
        
        kB = 8.617e-5 # where should I put this...

        # we want any input paramaters in results() to override self.params()
        params = lambda k: results(k) if results.has_key(k) else self.params(k)
        
        # initialize on hyperplane
        r = results("ReactionCoordinate")
        T = results("Temperature")
        
        self.initialize_hyperplane(r,0.)
        # TODO: check order in public PAFI
        self.run_script("PreRun")
        self.initialize_hyperplane(r,T)
        n_atoms = self.get_natoms()

        # the PAFI fix
        gamma = params("Friction")
        overdamped = params("OverDamped")
        seed = self.params.randint()
        cmd = f"fix pafi all pafi __pafipath {T} {gamma} {seed} "
        cmd += f"overdamped {overdamped} com 1\nrun 0"
        self.run_commands(cmd)

        # pre minimize (optional)
        if params("PreMin"):
            min_steps = params("MinSteps")
            self.run_commands(f"""
                min_style fire
                minimize 0 0.0001 {min_steps} {min_steps}
            """)
        if params("PostMin"):
            change_x = -self.gather("x",1,3)
        
        # PreThermalize (optional)
        self.run_script("PreTherm")
        results.set("MinEnergy",self.get_energy())
        
        # establish temperature time average and thermalize
        pre_therm_hp = self.extract_fix("pafi",size=4)
        steps = params("ThermSteps")
        ave_steps = params("ThermWindow")
        f_T = "c_pe" if overdamped==1 else "c_thermo_temp"
        self.run_commands(f"""
            reset_timestep 0
            fix __ae all ave/time 1 {ave_steps} {steps} {f_T}
            run {steps}
        """)
        sampleT = self.extract_fix("__ae")
        
        if overdamped==1:
            sampleT = (sampleT-results("MinEnergy"))/1.5/n_atoms/kB
        results.set("preT",sampleT)
        self.run_commands(f"""
            unfix __ae
            run 0""")
        
        # main sampling run
        steps = params("SampleSteps")
        self.run_commands(f"""
            reset_timestep 0
            fix __ae all ave/time 1 {steps} {steps} {f_T}
            """)
        if params("PostDump"):
            self.run_commands(f"""
            fix pafiax all ave/atom 1 {steps} {steps} x y z
            """)
        
        results = self.constrained_average(results)
        
        # get final temperature
        sampleT = self.extract_fix("__ae")
        if overdamped==1:
            sampleT = (sampleT-results("MinEnergy"))/1.5/n_atoms/kB
        results.set("postT",sampleT)
        self.run_commands(f"""
            unfix __ae
            run 0""")
        
        # average positions
        if params("PostDump"):
            dx = self.gather("f_pafiax",type=1,count=3)
            dx[:,0] -= self.gather("f_ux",type=1,count=1).flatten()
            dx[:,1] -= self.gather("f_uy",type=1,count=1).flatten()
            dx[:,2] -= self.gather("f_uz",type=1,count=1).flatten()
            dx = self.pbc(dx)
            results.set("MaxDev",np.linalg.norm(dx,axis=1).max())
            results.set("Dev",dx.copy())
            del dx
            self.run_commands("unfix pafiax")
        
        # minimize to test for MaxJump
        if params("PostMin"):
            min_steps = params("MinSteps")
        else:
            min_steps = 1
        self.run_commands(f"""
            min_style fire
            minimize 0 0.0001 {min_steps} {min_steps}
        """)
        change_x += self.gather("x",1,3)
        results.set("MaxJump",self.pbc_dist(change_x,axis=1).max())
        results.set("Valid",bool(results("MaxJump")<params("MaxJumpThresh")))
        
        # unfix hyperplane
        self.run_commands("unfix pafi")
        self.run_script("PostRun")
        # rescale back.... not sure if this is required
        self.initialize_hyperplane(r,0.0)

        return results
  
            



        





        


