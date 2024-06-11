#!/bin/bash

# INPUT ANGLES IN DEGREES
THETA=74.1
PHI=35.0

for a in {0..4000..100}; do
    A=$( bc -l <<< $a*0.0001 )
    sbatch submit.polariton ${A} ${THETA} ${PHI}
done
