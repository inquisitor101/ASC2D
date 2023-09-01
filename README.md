# ASC2D
ASC2D stands for `a simple Cartesian 2D` solver. More specifically, this is a two-dimensional nodal discontinuous Galerkin (DG-FEM) solver. Its purpose is for testing different non-reflective boundary conditions (NRBCs) in the context of the non-linear Euler equations (EEs). Additional information and/or templates would be provided, if users are interested.

# Some Useful Information
* For documentation, please refer to [the documentation provided](https://github.com/inquisitor101/ASC2D/blob/main/doc.pdf).

* An example folder with pre-defined test-cases is found [here](https://github.com/inquisitor101/ASC2D/tree/main/example).

* Additional tools for post-processing are found [here](https://github.com/inquisitor101/ASC2D/tree/main/tools).

# Some Relevant Notes
* This has been tested on Ubuntu 20.04 and 22.04.

* The [GNU gcc](https://gcc.gnu.org/) compiler has to be version 8 or higher.

* In order to compile the code, the [Make](https://www.gnu.org/software/make/) build utility is required.

* During compilation and running, the [Eigen](https://eigen.tuxfamily.org) library is needed.

* In order to run in parallel (via shared memory), an [OpenMP](https://www.openmp.org/) library must be available.

* For visualizing the solution, the open-source visualization software [ParaView](https://www.paraview.org/) is recommended.

# In a Nutshell
* To download, clone this repository via: `git clone  https://github.com/inquisitor101/ASC2D`
* To set-up Eigen library submodule, run: `git submodule init` followed by: `git submodule update`
* To install, enter `asc/` and run: `make`
* To run any test case, make sure to create the below eight directories:

 `data/` `anim/`  `proc/`  `logs/` `gnuplot/` `surf/` `zone/` `init/`
* To specify the number of threads/cores (`nThreads`), use: `export OMP_NUM_THREADS=nThreads`
* To run the solver, use: `/path/to/executable/ACS3 paramDict.cfg > logs/info.log 2>&1 &`
* To clean the simulation directory, simply execute: `./sweeper.clean`
* ... and most importantly: **HAVE FUN!!!**



