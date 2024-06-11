#!/bin/bash


for state in {0..100}; do
    /scratch/bweight/software/bader_charge/bader dens.${state}.cube
    mv AVF.dat AVF_${state}-${state}.dat
    mv ACF.dat ACF_${state}-${state}.dat
    mv BCF.dat BCF_${state}-${state}.dat
done