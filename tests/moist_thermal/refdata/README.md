stats_ens_1000.txt is a file with mean, std dev, min and max calculated from 1000 moist thermal simulations.
Differences between runs are caused by stochasticity in Lagrangian microphysics.
Eulerian microphysics (1-mom bulk) have no stochasticity, so all simulations give the same results.

However, Eulerian microphysics give slightly different results with MPI than without MPI (why?).
Therefore a file with mean, std dev, min and max calculated from 2 runs, 1 with MPI and 1 without MPI, is created (stats_ens_2.txt).

It is not known if Lagrangian microphysics give different results with MPI and without MPI.
With MPI, the results are within the spread of simulations without MPI.

Lagrangian microphyscis are compared with stats_ens_1000.txt.
Eulerian microphyscis are compared with stats_ens_2.txt.

