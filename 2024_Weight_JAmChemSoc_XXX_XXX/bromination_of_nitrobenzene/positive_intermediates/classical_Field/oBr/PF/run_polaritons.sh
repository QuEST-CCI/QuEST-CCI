#!/bin/bash

# INPUT ANGLES IN DEGREES
THETA=90.0
PHI=90.0
WC=2.8
WCAU=$( bc -l <<< ${WC}/27.2114 )

for a in {-150..150..5}; do
    A=$( bc -l <<< $a/10/514/$WCAU )
    sbatch submit.polariton ${A} ${WC} ${THETA} ${PHI}
done
