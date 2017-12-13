#!/bin/sh
totstep=720
nstep=30
step=0
while [ $step -le $totstep ]
do
  step_next=$(($step+$nstep))
# + $nstep
  echo $step, $step_next
  /home/piotr/praca/ParaView-5.4.1-Qt5-OpenGL2-MPI-Linux-64bit/bin/pvbatch --use-offscreen-rendering /home/piotr/praca/paraview_scripts/DYCOMS_animation.py $step $step_next
  step=$(($step_next+1))
done
