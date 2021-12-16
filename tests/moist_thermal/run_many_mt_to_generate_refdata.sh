# to generate moist thermal refdata 
# run the test multiple times using this script (from build/moist_thermal directory)
# then use moist_thermal_create_refdata to get mean and std_dev of statistics from these runs

for i in {1..40}
do
  mkdir -p $i
  cd $i
  OMP_NUM_THREADS=1 ../moist_thermal ../.. &
  cd ..
  sleep 1s # sleep to get different rng seed in next run
done
