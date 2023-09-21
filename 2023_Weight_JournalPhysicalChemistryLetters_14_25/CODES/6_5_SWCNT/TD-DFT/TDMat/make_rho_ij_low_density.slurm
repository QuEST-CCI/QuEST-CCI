#!/bin/bash
#SBATCH -p action
#SBATCH -J TransDens
#SBATCH -o output.slurm
#SBATCH -t 5-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem 10GB

cp ~/THUNDER_HOME/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn_3.6_bin_Linux_noGUI/settings.ini .

for state in {19..20}; do

echo "Generating first transition density for transition state ${state}."

~/THUNDER_HOME/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn << EOF
../geometry.fchk
18
1
../geometry.out
${state}
1
1
13
0
0
0
0
0
0
0
EOF

mv transdens.cub trans-${state}.cube

done
