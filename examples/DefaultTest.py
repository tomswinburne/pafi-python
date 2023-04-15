import sys,glob
from mpi4py import MPI
sys.path.insert(1,'../')
from pafi import PAFIParser,PAFIManager,ResultsProcessor

# default parameter test, loading from directors
params = PAFIParser()

# find and sort file list, assign therm to PAFI parameters
fl = glob.glob("systems/EAM-EAM-SIA-Fe/image_*")
params.set_path(sorted(fl,key=lambda x:x.split("_")[-1]))

# Set potential
params.set_potential("systems/EAM-EAM-SIA-Fe/Fe.eam.fs")


# Print all parameter information
print(params.info())

if MPI.COMM_WORLD.Get_rank()==0:
    print(params.info())

manager = PAFIManager(MPI.COMM_WORLD,params=params)
manager.run()
manager.close()

exit()



