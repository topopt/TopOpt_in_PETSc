TopOpt_in_PETSc
===============
A 3D large-scale topology optimization code using PETSc
===============

The code (or framework) presented on this page is a fully parallel framework for conducting very large scale topology optimziation on structured grids. For more
details see www.topopt.dtu.dk/PETSc.

Updated and refactored to remove dependence of TopOpt.cc/h in all other classe,
June, 2019, Niels Aage

To clone repository:
>> git clone https://github.com/topopt/TopOpt_in_PETSc.git

NOTE: The code requires PETSc version 3.11.0 or newer ! Also note that the code is not tested against the development branch on git.

This code has been tested on:
- Linux systems including: Ubuntu 18.04, Red hat enterprise linux 8

This code requires the following external software to work:
- PETSc version 3.11.4 or earlier (though never than 3.8.x)
- Requires LAPACK/BLAS
- Requires MPI

Compile following rules in makefile_ref

Normal compilation time of framework, e.g. 4s: "make topopt -j"

Run the base example by typing e.g.: "mpirun -np 4 ./topopt"

Postprocess results using Python 2.6: "bin2vtu #" where # refers to the iteration number

Visualize using ParaView (version 5.7 or earlier)

The expected result of the base code is the (but on a coarse mesh!) cantilever beam from:
Aage, N., Andreassen, E., & Lazarov, B. S. (2015). Topology optimization using PETSc: An easy-to-use, fully parallel, open source topology optimization framework. Structural and Multidisciplinary Optimization, 51(3), 565–572. https://doi.org/10.1007/s00158-014-1157-0

Extensions: 
An extension of the code including manufacturing filters/constraints can be found here:
https://github.com/edofersan/MaximumSize_on_TopOpt_in_PETSc



