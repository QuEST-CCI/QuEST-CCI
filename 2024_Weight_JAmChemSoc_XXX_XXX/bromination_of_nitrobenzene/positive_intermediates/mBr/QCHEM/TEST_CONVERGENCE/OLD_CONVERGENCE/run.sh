#!/bin/bash

A0=0.4

NF=5
for NM in {2,3,4,5,10,25,50,75,100,200,300,400,500,600,700,800,900,1000}; do
    echo $NM $NF
    #rm -r PF_NF${NF}_NM${NM}/
    mkdir -p PF_NF${NF}_NM${NM}/
    cd PF_NF${NF}_NM${NM}/
        cp ../Pauli-Fierz_DdotE.py .
        cp ../submit.polariton .
        sed -i "s/BASH_NM/$NM/g" Pauli-Fierz_DdotE.py
        sed -i "s/BASH_NF/$NF/g" Pauli-Fierz_DdotE.py
        sbatch submit.polariton ${A0} 5.0
        cd ../
done

NM=50
for NF in {2,3,4,5,6,7,8,9,10,20,50}; do
    echo $NM $NF
    #rm -r PF_NF${NF}_NM${NM}/
    mkdir -p PF_NF${NF}_NM${NM}/
    cd PF_NF${NF}_NM${NM}/
        cp ../Pauli-Fierz_DdotE.py .
        cp ../submit.polariton .
        sed -i "s/BASH_NM/$NM/g" Pauli-Fierz_DdotE.py
        sed -i "s/BASH_NF/$NF/g" Pauli-Fierz_DdotE.py
        sbatch submit.polariton ${A0} 5.0
        cd ../
done


