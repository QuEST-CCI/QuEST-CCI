#!/bin/bash

# INPUT ANGLES IN DEGREES
THETA=69.4 # o/m 74.1, p/m 69.4
PHI=79.9 # o/m 35.0, p/m 79.9

for a in {0..4000..100}; do
    A=$( bc -l <<< $a*0.0001 )
    sbatch submit.polariton ${A} ${THETA} ${PHI}
done
