import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

DIP_TEMP = np.loadtxt("dipole_matrix_E.dat")

MU = np.zeros((21,21,3))

for line in DIP_TEMP:
    j = int( line[0] )
    k = int( line[1] )
    MU[j,k,:] = np.array( line[2:], dtype=float )
    MU[k,j,:] = np.array( line[2:], dtype=float )

np.savetxt( "MU_z.dat", np.abs(MU[:,:,2]) )

#cmap = mpl.cm.hot
#cmap = mpl.cm.hot_r
#cmap = mpl.cm.Greys
#cmap = mpl.cm.tab20c
#cmap = mpl.cm.rainbow
#cmap = mpl.cm.terrain
cmap = mpl.cm.terrain_r
plt.imshow( np.abs(MU[:,:,-1]) , origin='lower', cmap=cmap, norm=mpl.colors.Normalize(vmin=0, vmax=30))

plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=30), cmap=cmap))

#plt.xlim(0,21)
#plt.ylim(0,21)
plt.xlabel("State Index", fontsize=10)
plt.ylabel("State Index", fontsize=10)
#plt.title(f"Absorption (NExc:{Nad} Nfock: {NFock})",fontsize=10)
plt.tight_layout()
plt.savefig(f"Dipole_Matrix_Ez.jpg", dpi=600)
plt.clf()
