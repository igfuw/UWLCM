To generate moist thermal refdata:

1. run the test multiple times using run_many_mt_to_generate_refdata.sh (from build/moist_thermal directory)

2. use moist_thermal_create_refdata (e.g. ./moist_thermal_create_refdata {1..100}) to get mean and std_dev of statistics from these runs

3. copy files with std_Dev and mean to refdata dir
