Example code for dynamical mean-field theory (DMFT) calculation using the NRG Ljubljana code as
the impurity solver

Model: Kondo lattice model (KLM), S=1/2; model description in template/.

Features:
- use of templates (i.e., Mathematica is not required at run time)
- improved NRG discretization scheme (Zitko, Pruschke, 2008)
- improved estimator for the self-energy (Kugel, 2022)
- support for arbitrary density of states (tabulated in file DOS.dat)
- robust band occupancy control by shifting chemical potential mu
- simple mixing (implemented at the level of hybridisation function)
- adaptive grid for better capturing sharp spectral features
- transport calculation using external "bubble" code

Requirements:
- NRG Ljubljana with associated tools (hilb, kk, adapt, nrgchain, broaden, resample, matrix, diag, unitary)
- associated scripts (getparam, scaley, getiter, newiter, subtracty...), in github repo rokzitko/nrgljubljana under scripts/.
- perl
- Python (for occupancy_control)
- m4 macro processor
- bubble (optional)

Two modes of operation:
- local: a script named "mynrgrun" must exist; a minimal version just calls "nrg", but typically you will want to set up
  the environment (number of threads, working directory, piping of output to files); an example is provided
- slurm: a script named "subslurm" must exist to create a job script and submit it to the cluster for execution

Rok Zitko, 2026
