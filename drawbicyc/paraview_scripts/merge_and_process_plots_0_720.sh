#!/bin/bash
totstep=720
nstep=30
step=0
iter=0
while [ $step -le $totstep ]
do
  step_next=$(($step+$nstep))
# + $nstep
  echo $iter, $step, $step_next
  arr[$iter]=/home/piotr/praca/wyniki/DYCOMS3D/a02_SD30_dt1_4_12_17_outfreq30/paraview_batch_animation_${step}_${step_next}.avi
  #arr[0]=0
  step=$(($step_next+1))
  ((iter++))
done

rm /home/piotr/praca/wyniki/DYCOMS3D/a02_SD30_dt1_4_12_17_outfreq30/paraview_batch_animation_merged_0_720.avi
avimerge -o /home/piotr/praca/wyniki/DYCOMS3D/a02_SD30_dt1_4_12_17_outfreq30/paraview_batch_animation_merged_0_720.avi -i ${arr[*]} 

#convert the merged file into an mp4
rm /home/piotr/praca/wyniki/DYCOMS3D/a02_SD30_dt1_4_12_17_outfreq30/paraview_batch_animation_merged_0_720.mp4
avconv -i /home/piotr/praca/wyniki/DYCOMS3D/a02_SD30_dt1_4_12_17_outfreq30/paraview_batch_animation_merged_0_720.avi -c:v libx264 -c:a copy /home/piotr/praca/wyniki/DYCOMS3D/a02_SD30_dt1_4_12_17_outfreq30/paraview_batch_animation_merged_0_720.mp4

# change frame rate of the mp4
#ffmpeg -i /home/piotr/praca/wyniki/DYCOMS3D/a02_SD30_dt1_4_12_17_outfreq30/paraview_batch_animation_merged_0_720.mp4 -r 16 -filter:v "setpts=0.25*PTS" /home/piotr/praca/wyniki/DYCOMS3D/a02_SD30_dt1_4_12_17_outfreq30/paraview_batch_animation_merged_0_720.mp4

# use interpolation to increase frame rate (smooth it out)
rm /home/piotr/praca/wyniki/DYCOMS3D/a02_SD30_dt1_4_12_17_outfreq30/paraview_batch_animation_merged_0_720_smoothto60.mp4
ffmpeg -i /home/piotr/praca/wyniki/DYCOMS3D/a02_SD30_dt1_4_12_17_outfreq30/paraview_batch_animation_merged_0_720.mp4 -filter "minterpolate='mi_mode=mci:mc_mode=aobmc:vsbmc=1:scd=none:fps=60'" /home/piotr/praca/wyniki/DYCOMS3D/a02_SD30_dt1_4_12_17_outfreq30/paraview_batch_animation_merged_0_720_smoothto60.mp4
