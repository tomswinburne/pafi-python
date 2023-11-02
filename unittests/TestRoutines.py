import sys
sys.path.insert(1,'../')
from mpi4py import MPI

def test_plotting():
    from pafi import ResultsProcessor
    p = ResultsProcessor(data_path="dumps/pafi_data_0.csv")
    print(p.integrate())

def test_complete_input():
    """Load in complete configuration file
    """
    from pafi import PAFIManager,ResultsProcessor
    manager = PAFIManager(MPI.COMM_WORLD,"./CompleteConfiguration.xml")
    manager.run()
    manager.close()
    p = ResultsProcessor(data_path=manager.parameters.csv_file)
    print(p.integrate())


def test_partial_input():
    """Overwrite default parameters with partial configuration file
    """
    from pafi import PAFIManager
    manager = PAFIManager(MPI.COMM_WORLD,"./PartialConfiguration.xml")
    manager.run()
    manager.close()

def test_python_input():
    """Overwrite parameters within python script
    """
    from pafi import PAFIManager,PAFIParser
    # load in default parameters
    parameters = PAFIParser()
    # set path wildcard, assuming some_file_name_INTEGER.dat format
    parameters.set_pathway("systems/EAM-EAM-Vac-Fe/image_*.dat")
    # set interatomic potential
    parameters.set_potential("systems/EAM-EAM-Vac-Fe/Fe.eam.fs")
    # restrict to zero temperature
    parameters.axes["Temperature"] = [100.]
    parameters.set("nRepeats",2)
    parameters.set("SampleSteps",10)
    parameters.set("ThermSteps",10)
    parameters.set("ThermWindow",10)
    
    manager = PAFIManager(MPI.COMM_WORLD,parameters=parameters)#,restart_data)
    manager.run()
    manager.close()

def test_custom_worker():
    """Custom worker for thermodynamic integration
    """
    from pafi import PAFIManager,PAFIParser,PAFIWorker,ResultsHolder
    from typing import List
    class CustomPAFIWorker(PAFIWorker):
        def __init__(self, comm: MPI.Intracomm, 
                     parameters: PAFIParser, 
                     tag: int, rank: int, roots: List[int]) -> None:
            super().__init__(comm, parameters, tag, rank, roots)

        def constrained_average(self,results:ResultsHolder)->ResultsHolder:
            """Perform constrained PAFI average for thermodynamic integration
            In CustomConfiguration.xml we establish a new axis, 'Lambda', 
            and have a set of scripts to set the potential as 
            V = Lambda*(V_1-V_2) + V_2
            then evaluate <V_1-V_2 | Lambda > 
            """
            parameters = lambda k: results(k)\
                if results.has_key(k) else self.parameters(k)
    
            steps = parameters("SampleSteps")
        
            fixname = self.setup_pafi_average(steps,"avepafi")

            # set average fix, run for SampleSteps
            self.run_commands(f"""
            fix dV all ave/time 1 {steps} {steps} v_dV
            run {steps}
            """)

            # extract average
            ave_dV = self.extract_fix("dV",size=1)
            results.set("ave_dV",ave_dV)

            # remove fix
            self.run_commands("unfix dV")
            
            results = self.extract_pafi_data(results,fixname)
            return results

    manager = PAFIManager(MPI.COMM_WORLD,
                          xml_path="./CustomConfiguration.xml",Worker=CustomPAFIWorker)
    manager.run(print_fields=\
                ["Temperature","Lambda","ReactionCoordinate","aveF","ave_dV"])

if __name__ == "__main__":
    #test_plotting()
    test_python_input()
    #test_partial_input()
    #test_custom_worker()

exit()



