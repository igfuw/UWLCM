stats_ens_1000.txt is a file with mean, std dev, min and max calculated from 1000 moist thermal simulations (without MPI).
Differences between runs are caused by stochasticity in Lagrangian microphysics.
Eulerian microphysics (1-mom bulk) have no stochasticity, so all simulations give the same results.

However, Eulerian microphysics give slightly different results with MPI than without MPI (why? This happens even for mpirun with 1 process (?)).
Therefore a separate file with min and max from one blk_1m MPI run is created (stats_mpi_blk_ens_1.txt) and used for comparison in blk_1m MPI.

It is not known if Lagrangian microphysics give different results with MPI and without MPI.
With MPI, the results are within the spread of simulations without MPI.

