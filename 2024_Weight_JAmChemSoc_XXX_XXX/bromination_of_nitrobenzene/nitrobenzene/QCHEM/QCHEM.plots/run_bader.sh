#!/bin/bash

for cube_file in $(find . -type f -name "*cube"); do
    /scratch/bweight/software/bader_charge/bader ${cube_file}
    mv AVF.dat ${cube_file}.AVF
    mv ACF.dat ${cube_file}.ACF
    mv BCF.dat ${cube_file}.BCF
done
