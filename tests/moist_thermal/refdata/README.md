stats_lgrngn_ens_1000.txt - mean, std dev, min and max calculated from 1000 moist thermal simulations with Lagrangian microphysics, without MPI.
Differences between runs are caused by stochasticity in Lagrangian microphysics.

Eulerian microphysics (1-mom bulk) have no stochasticity, so all simulations give the same results.
stats_blk_ens_1.txt - results of blk_1m runs without MPI
stats_mpi_blk_ens_1.txt - results of blk_1m runs with MPI
Why are there (small) differences between MPI and non-MPI runs?
This happens even for mpirun with 1 process (?).

It is not known if Lagrangian microphysics give different results with MPI and without MPI.
With MPI, the results are within the spread of simulations without MPI.

