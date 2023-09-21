import numpy as np
import subprocess as sp
import matplotlib as mpl
from matplotlib import pyplot as plt
import os

NExcStates = 500 # Number of excited states S1, S2, S3, ...

dipole_mat = np.zeros(( NExcStates+1, NExcStates+1, 3 ))

# Read in transition dipole moments
transDipFile = open("transdipmom.txt","r").readlines()
for line in transDipFile:
    t = line.split()
    if (  len(t) == 7 and t[0] != "i" and t[0] != "Transition" ):
        i = int( t[0] )
        j = int( t[1] )
        if ( i <= NExcStates and j <= NExcStates and i != j ): # Exclude terms from ground state and permanant dipoles
            dxyz = np.array( t[2:5] )
            dipole_mat[ i, j ] = dxyz
            dipole_mat[ j, i ] = dxyz

NPlot = 20
cmap = mpl.cm.hot
for d in [0,1,2]:
    plt.imshow( np.abs(dipole_mat[:NPlot,:NPlot,d]) , origin='lower', cmap=cmap, norm=mpl.colors.Normalize(vmin=0, vmax=3.0))
    plt.colorbar()
    #plt.title("µ (a.u.)")
    plt.tight_layout()
    if ( d == 0 ):
        plt.savefig(f"Dipole_Matrix_Ex_{NPlot}.jpg", dpi=600)
    if ( d == 1 ):
        plt.savefig(f"Dipole_Matrix_Ey_{NPlot}.jpg", dpi=600)
    if ( d == 2 ):
        plt.savefig(f"Dipole_Matrix_Ez_{NPlot}.jpg", dpi=600)
    plt.clf()








