import numpy as np
import subprocess as sp
import os
from glob import glob
from matplotlib import pyplot as plt

CHI = np.arange(-15,15,0.5) # V/nm --> a.u.

NCHI  = len( CHI )
ENERGY  = np.zeros( (NCHI,3) ) #(x,y,z)

for dind,d in enumerate(["X","Y","Z"]):
    os.chdir(f"DATA/")
    for step in range( NCHI ):
        print( f"\tStep {step+1} of {NCHI}" )
        os.chdir(f"STEP_{step}/")
        ENERGY[step,dind] = float( sp.check_output("grep 'SCF Done' geometry.out | tail -n 1 | awk '{print $5}'",shell=True) )
        ENERGY[step,dind] *= 630
        os.chdir(f"../")
    os.chdir(f"../")

np.savetxt( "Energy.dat", np.c_[CHI, ENERGY[:,0], ENERGY[:,1], ENERGY[:,2]], fmt="%1.5f")

plt.plot( CHI, ENERGY - ENERGY[NCHI//2] )
plt.xlabel("Field Strength, $\mathcal{E}$ (V/nm)", fontsize=15)
plt.ylabel("Energy, (kcal/mol)", fontsize=15)
plt.savefig("Energy.jpg", dpi=300)

