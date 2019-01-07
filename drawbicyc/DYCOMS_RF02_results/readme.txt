GCSS Boundary Layer Cloud Working Group Results for DYCOMS-II RF02 
Intercomparison

All the final model output used in the intercomparison are provided in 
netCDF format.  Additional processing not applied to the provided output 
data was performed for the plots and analysis as follows. The scalar time 
series data were linearly interpolated to a temporal grid with uniform 
spacing of 300 seconds. The profiles were linearly interpolated to a 
vertical grid with uniform spacing of 2 meters, and where profile 
variables were missing but needed at the surface they were linearly 
extrapolated from the overlying two layers.  Profiles were then scaled by 
the corresponding half-hour average inversion height (from the scalar time 
series) before being temporally averaged.

The array named "missing" is set to zero unless a particular model 
configuration is missing, in which case it is set to one (e.g., for 
COAMPS_SL without drizzle but with cloud-water sedimentation).

The dimension named "drizzle" has a length of 2, with the first (second) 
element corresponding to drizzle being off (on).  Likewise, the first 
(second) element of the dimension named "sed" corresponds to cloud-water 
sedimentation being omitted (included).

The character array named "zmap" gives the name of each variable and the 
name of the vertical array the variable is located (e.g., the pair "ql" 
and "zt" for a simulation means that the ql variable for that simulation 
is located at the zt grid points).

NetCDF fill values are used for all missing data.

The provided Fortran-90 demonstration program reads some time series and 
profile data, then computes and prints out liquid water path averaged over 
the last four hours from the time series and from the profiles.  Note that 
for 6 of the 14 simulations (with drizzle and cloud-water sedimentation 
included) the numbers do not match because of inconsistencies in the 
submitted data.  For the analysis in the manuscript, the maximum of the 
two was used where they disagree.
