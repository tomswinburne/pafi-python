"""

    An example of how we can customize the PAFI run "in place"

    Clearly, it is better (for readability) to have the scripts
    in a separate XML file!

"""

import numpy as np
import os,sys
import pandas as pd
from mpi4py import MPI
sys.path.insert(1,'../')
from pafi import PAFIManager,PAFIWorker,Parser,ResultsHolder

params = Parser("./StandardConfiguration.xml")
params.parameters["LogLammps"] = True
params.scripts["Input"] = """
    units metal
    atom_style atomic
    atom_modify map array sort 0 0.0
    neigh_modify every 2 delay 10 check yes page 1000000 one 100000
    read_data  ../examples/systems/EAM-EAM-SIA-Fe/image_1.dat
    pair_style    hybrid/overlay eam/fs eam/fs
    pair_coeff * * eam/fs 1 ../examples/systems/EAM-EAM-SIA-Fe/Fe.eam.fs Fe
    pair_coeff * * eam/fs 2 ../examples/systems/EAM-EAM-SIA-Fe/Fe_mm.eam.fs Fe
    run 0
"""
params.scripts["PreRun"] = """
    variable l1 equal 1.0-%Lambda%
    variable l2 equal 0.0+%Lambda%
    run 0
    fix escale1 all adapt 1 pair eam/fs:1 scale * * v_l1
    fix escale2 all adapt 1 pair eam/fs:2 scale * * v_l2
    run 0
    compute V_1 all pair eam/fs 1
    compute V_2 all pair eam/fs 2
    variable V_1 equal (c_V_1/(0.000001+v_l1))
    variable V_2 equal (c_V_2/(0.000001+v_l2))
    variable dV equal v_V_2-v_V_1
"""
params.scripts["PostRun"] = """
    variable dV delete
    uncompute V_1
    uncompute V_2
    unfix escale1
    unfix escale2
"""

params.axes["Lambda"] = np.linspace(0.,1.0,5)


class CustomPAFIWorker(PAFIWorker):
    def __init__(self, comm: MPI.Intracomm, params: Parser, tag: int) -> None:
        super().__init__(comm, params, tag)
    
    def constrained_average(self,results:ResultsHolder)->ResultsHolder:
        
        # ensure results inputs should override self.params()
        params = lambda k: results(k) if results.has_key(k) else self.params(k)
        
        steps = params("SampleSteps")
        
        fixname = self.setup_pafi_average(steps,"avepafi")
        
        self.run_commands(f"""
            fix dV all ave/time 1 {steps} {steps} v_dV
            run {steps}
        """)
        ave_dV = self.extract_fix("dV",size=1)
        results.set("ave_dV",ave_dV)

        results = self.extract_pafi_data(results,fixname)
        return results

manager = PAFIManager(MPI.COMM_WORLD,params=params,Worker=CustomPAFIWorker)
manager.run(print_fields=["Rank","Temperature","Lambda","ReactionCoordinate","aveF","ave_dV"])
manager.close()
exit()



