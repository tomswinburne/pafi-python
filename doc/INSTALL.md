![](pafi_title.png)
## [Main Page](README.md)

# Installation

## Compile `LAMMPS` with `EXTRA-FIX` package

1. `PAFI` is now integrated into `LAMMPS` as part of the `EXTRA-FIX` package, introduced in the 28 July 2021 Patch release. For earlier versions, use the `USER-MISC` package instead.
You can [download](https://lammps.sandia.gov/download.html) a tarball from the `LAMMPS`
website or clone the public repository with
```bash
git clone https://github.com/lammps/lammps.git
```

2. Install `EXTRA-FIX` and any packages you desire (e.g. replica for `NEB`)
```bash
cd /path/to/lammps/src
make yes-extra-fix ## yes-user-misc for versions before 28 July 2021 Patch release
make yes-replica # for NEB calculation
make yes-package # (e.g. manybody for EAM potentials etc)
```

3. In the appropriate Makefile add `-std=c++11` to `CCFLAGS` and `LINKFLAGS` and
add `-DLAMMPS_EXCEPTIONS` to `LMP_INC` to allow `PAFI` to read `LAMMPS` error messages.
This is very useful when running your own simulations. For `src/MAKE/Makefile.mpi` this reads
 ```make
CCFLAGS =	-g -O3 -std=c++11
LINKFLAGS =	-g -O3 -std=c++11
LMP_INC =	-DLAMMPS_GZIP -DLAMMPS_MEMALIGN=64  -DLAMMPS_EXCEPTIONS
```

4. Compile static library and binary Consult [LAMMPS documentation](http://lammps.sandia.gov/doc/Section_start.html) for details
```bash
   make mpi mode=lib # liblammps_mpi.a library for pafi
   make mpi # lmp_mpi binary for running initial NEB calculation if desired
```

4. Copy library to your local lib/ and headers to local include/, at e.g. ${HOME}/.local
```bash
  export PREFIX=${HOME}/.local # example value
  cp liblammps_mpi.a ${PREFIX}/lib/liblammps_mpi.a
  mkdir ${PREFIX}/include/lammps
  cp *.h ${PREFIX}/include/lammps/
```

5. You should be able to run `mpirun -np 4 python test.py` without errors, where `test.py` reads
```python
from mpi4py import MPI
from lammps import lammps
lmp = lammps(comm=MPI.COMM_WORLD)
lmp.close()
```

## Install `PAFI`
Happily, `PAFI` is pure python and so does not require compilation. 
You can use `PAFI` in your code simply with three lines:
```python
import sys
sys.path.insert(1,/path/to/pafi/root)
import pafi
```
Alternatively, the `PAFI` library can be installed with `pip`:
```bash
cd /path/to/pafi/root
pip install -e .
```


