import numpy as np
from matplotlib import pyplot as plt

wc_list = np.arange( 4.5, 10.0, 0.05 )

EDiff = np.zeros(( len(wc_list), 2, 2 ))
for wcIND, wc in enumerate( wc_list ):
    EDiff[wcIND,:,:] = np.loadtxt(f"image_data/EPol_Diff__Ez_A00.04_wc{round(wc,5)}.dat")

AC_Diff = np.zeros(( len(wc_list) ))
R_AC = np.zeros(( len(wc_list) ))
for wcIND in range( len(wc_list) ):
    if ( EDiff[wcIND,0,1] < EDiff[wcIND,1,1] ):
        R_AC[wcIND]    = EDiff[wcIND,0,0]
        AC_Diff[wcIND] = EDiff[wcIND,0,1]
    else:
        R_AC[wcIND]    = EDiff[wcIND,1,0]
        AC_Diff[wcIND] = EDiff[wcIND,1,1]

print( np.c_[R_AC, AC_Diff] )

#plt.plot( wc_list, EDiff[:,0,1], "-" )
#plt.plot( wc_list, EDiff[:,1,1], "-" )
plt.plot( wc_list, AC_Diff, "o-", label="A$_0$ = 0.04" )

plt.legend()

plt.xlabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.ylabel("Avoided Crossing Magnnitude (meV)",fontsize=15)

plt.xlim(4.5,10)
plt.ylim(0)

plt.tight_layout()
plt.savefig("Avoided_Crossing_vs_wc.jpg")
np.savetxt( "Avoided_Crossing_vs_wc.dat", np.c_[wc_list, AC_Diff] )

