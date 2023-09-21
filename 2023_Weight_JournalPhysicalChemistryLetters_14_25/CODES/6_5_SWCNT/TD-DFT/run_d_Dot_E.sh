#!/bin/bash

for i in {0..250..1}; do # Loop from 0.005 to 0.250
	for j in {100..200}; do # Loop from 8.00 to 8.50
		eta=$( bc -l <<< $i*0.001 )
		wc=$( bc -l <<< $j*0.01 )
		echo;
		echo "eta = ${eta}"
		echo "wc = ${wc} eV"
        	python3 ../../../d_dot_E_Braden.py ${eta} ${wc} 
        	#python3 plot_polaritons.py ${eta} ${wc} 
	done
done

