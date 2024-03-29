import sys,glob
import numpy as np
sys.path.insert(1,'../')
from mpi4py import MPI
import argparse

def test_integration():
    from pafi import ResultsProcessor
    """
        All test data here is for the same system, currently at zero temperature
        TODO: have finite temperature tests..
    """
    csv_list = glob.glob("dumps/pafi_data_*.csv")
    p = ResultsProcessor(data_path=csv_list)
    
    # Returns a pandas DataFrame and a list of dictionaries
    x_key = 'ReactionCoordinate'
    y_key = 'FreeEnergyGradient'
    _ , integrated_data = p.integrate(argument=x_key,
                    target=y_key,remesh=10,
                    return_remeshed_array=True)
    for int_dat in integrated_data:
        ii = int_dat[y_key+"_integrated"].argmax()
        f_b = np.round(int_dat[y_key+"_integrated"][ii],4)
        f_e = np.round(int_dat[y_key+"_integrated_err"][ii],4)    
        print(f"""
            Temperature: {int_dat['Temperature']}
            Free Energy Barrier: {f_b} +- {f_e}eV
            """)
    
        
def test_complete_input():
    """Load in complete configuration file
    """
    from pafi import PAFIManager
    config = "./configuration_files/CompleteConfiguration_TEST.xml"
    manager = PAFIManager(MPI.COMM_WORLD,config)
    manager.run()
    manager.close()
    
def test_partial_input():
    """Overwrite default parameters with partial configuration file
    """
    from pafi import PAFIManager
    config = "./configuration_files/PartialConfiguration_TEST.xml"
    manager = PAFIManager(MPI.COMM_WORLD,config)
    manager.run()
    manager.close()

def test_python_input():
    """Overwrite parameters within python script
    """
    from pafi import PAFIManager,PAFIParser
    # load in default parameters
    parameters = PAFIParser()
    # set path wildcard, assuming some_file_name_INTEGER.dat format
    parameters.set_pathway("systems/EAM-SIA-Fe/image_*.dat")
    # set interatomic potential
    parameters.set_potential("systems/EAM-SIA-Fe/Fe.eam.fs")
    # restrict to zero temperature
    parameters.axes["Temperature"] = [100.]
    parameters.set("nRepeats",2)
    parameters.set("SampleSteps",10)
    parameters.set("ThermSteps",10)
    parameters.set("ThermWindow",10)
    
    manager = PAFIManager(MPI.COMM_WORLD,parameters=parameters)
    manager.run()
    manager.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
    PAFI test routines: 
        Usage: 
            # e.g. 4 cores for sampling
            mpirun -np 4 python TestRoutines.py -t complete
            mpirun -np 4 python TestRoutines.py -t partial
            mpirun -np 4 python TestRoutines.py -t python

            # just test postprocessing
            python TestRoutines.py -t integrate
            """)
    
    options =  ['complete','partial','python','integrate']
    
    parser.add_argument('-t', '--test', help='Must be in '+" ".join(options))
    args = parser.parse_args()
    test = args.test.strip() 
    assert test in options
    
    if test==options[0]:
        test_complete_input()
    elif test==options[1]:
        test_partial_input()
    elif test ==  options[2]:
        test_python_input()
    elif test==options[3]:
        test_integration()
    

exit()



