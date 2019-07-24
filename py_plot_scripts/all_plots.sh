#!/bin/bash

FOLDER="output_lgr/dycoms"
SUB_FOLDER="6hrs_highres/exp"

for NUM in {1..5}
do
        echo starting exp$NUM
        python3 plot_series.py $FOLDER $SUB_FOLDER$NUM
        python3 plot_profs.py $FOLDER $SUB_FOLDER$NUM
        python3 plot_size_dist_new.py $FOLDER $SUB_FOLDER$NUM
        python3 plot_2Dmap.py $FOLDER $SUB_FOLDER$NUM
        echo done with exp$NUM
done
