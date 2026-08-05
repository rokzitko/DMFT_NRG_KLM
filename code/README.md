How to start:

- Copy param.mu.init (initial value for chemical potential mu) to param.mu. Note that param.mu is not in the repo, because it is variable (occupancy control).
- If DOS.dat (density of states of the band) does not exist, create it by running mkDOS.
- Modify model and solver settings in param.loop if needed.
- Start calculation by running DMFT scripts in current directory.
