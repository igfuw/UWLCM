#!/bin/bash

FOLDER="../output"
SUB_FOLDER="mu_0.01e-6_kappa_1.3/exp"

for NUM in {1..5}
do
        echo starting exp$NUM
        python3 plot_series.py $FOLDER $SUB_FOLDER$NUM
        python3 plot_profs.py $FOLDER $SUB_FOLDER$NUM
        python3 plot_size_dist.py $FOLDER $SUB_FOLDER$NUM
        #python3 plot_2Dmap.py $FOLDER $SUB_FOLDER$NUM
        echo done with exp$NUM
done
