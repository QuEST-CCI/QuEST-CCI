import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

DIP_TEMP = np.loadtxt("dipole_matrix_E.dat")

MU = np.zeros((21,21,3))

for line in DIP_TEMP:
    j = int( line[0] )
    k = int( line[1] )
    if ( j >= len(MU)-1 or k >= len(MU)-1 ):
        continue
    MU[j,k,:] = np.array( line[2:], dtype=float )
    MU[k,j,:] = MU[j,k,:] * 1.0

cmap = mpl.cm.hot

np.savetxt("MU_x.dat", MU[:,:,0])
np.savetxt("MU_x_ABS.dat", np.abs( MU[:,:,0] ))
np.savetxt("MU_norm_ABS.dat", np.abs( np.linalg.norm(MU[:,:,:],axis=-1)))


#plt.imshow( np.abs(MU[:6,:6,0]) , origin='lower', cmap=cmap, norm=mpl.colors.Normalize(vmin=0, vmax=30))
plt.imshow( np.abs() , origin='lower', cmap=cmap, norm=mpl.colors.Normalize(vmin=0, vmax=30))

plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=30), cmap=cmap), pad=0.01, )

#plt.xlim(0,21)
#plt.ylim(0,21)
plt.xlabel("State Index", fontsize=10)
plt.ylabel("State Index", fontsize=10)
#plt.title(f"Absorption (NExc:{Nad} Nfock: {NFock})",fontsize=10)
plt.tight_layout()
plt.savefig(f"Dipole_Matrix_Ex.jpg", dpi=600)
plt.clf()
