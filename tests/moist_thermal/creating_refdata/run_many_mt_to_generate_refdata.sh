for i in {1..1000}
do
  mkdir -p $i
  cd $i
  OMP_NUM_THREADS=1 ../moist_thermal ../.. &
  cd ..
  sleep 1s # sleep to get different rng seed in next run
done
