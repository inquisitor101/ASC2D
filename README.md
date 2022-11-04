# ASC2D
ASC2D stands for `a simple Cartesian 2D` solver. More specifically, this is a two-dimensional nodal discontinuous Galerkin (DG-FEM) solver. Its purpose is for testing different non-reflective boundary conditions (NRBCs) in the context of the non-linear Euler equations. Additional information and/or templates would be provided, if users are interested.

# Some Useful Information
* For documentation, please refer to [the documentation provided](./doc.pdf) 
* An example folder with pre-defined test-cases is found [here](./asc/example/)

# In a Nutshell
* To install, run: `make init_git_submodules && make`
* To run any test case, make sure to create the directories: `data/`  `anim/`  `proc/`  `logs/`
* To specify the number of threads/cores, use: `export OMP_NUM_THREADS=nThreads`
* To run the solver, use: `/path/to/executable/ACS3 paramDict.cfg > logs/info.log 2>&1 &`
* To clean the simulation directory, simply execute: `./sweeper.clean`
* ... and most importantly: **HAVE FUN!!!**



