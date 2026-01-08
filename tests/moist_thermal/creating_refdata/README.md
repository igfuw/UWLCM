To generate moist thermal refdata:

1. run the test multiple times using run_many_mt_to_generate_refdata.sh (from build/moist_thermal directory)
   This will run large number (1000) of repetitions of both lgrngn and bulk, where bulk only needs one. This can be changed in common.hpp

2. use moist_thermal_create_refdata (e.g. ./creating_refdata/moist_thermal_create_refdata {1..100}) to get mean and std_dev of statistics from these runs.
   This will create refdata for both lgngn and bulk assuming same ensemble size. Modify common.hpp to change this.

3. copy files with std_Dev and mean to refdata dir. May require change of filename to add info about type of microphysics.
