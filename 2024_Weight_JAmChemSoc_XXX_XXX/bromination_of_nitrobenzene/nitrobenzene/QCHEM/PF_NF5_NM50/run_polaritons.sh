#!/bin/bash

for j in {0..5000..100}; do
    eta=$( bc -l <<< $j*0.0001 )
    #echo $eta
    sbatch submit.polariton ${eta}
done
