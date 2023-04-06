import numpy as np
import os,sys
import pandas as pd
from mpi4py import MPI
sys.path.insert(1,'../')
from pafi import PAFIManager
# initial test
manager = PAFIManager(MPI.COMM_WORLD,"TestConfiguraions.xml")
manager.run()
manager.close()
exit()



