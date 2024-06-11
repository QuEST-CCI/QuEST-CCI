import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import os
from scipy.interpolate import interp2d
import subprocess as sp

NF  = 5
NM  = 50
A0_LIST = np.arange( 0.0, 0.4+0.02, 0.02 )
#A0_LIST = np.array([ 0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
WC_LIST = np.arange( 0.0, 4.1, 0.1 )
#WC_LIST = np.array([ 0.0, 5.0, 10.0, 15.0, 20.0])
DATA_DIR = "PLOTS_DATA_CRITICAL_ANGLE/"

THETA = 74.1 # 74.1/69.4 oBr/mBr
PHI   = 35.0 # 35.0/79.9 pBr/mBr

NPOL = NM*NF
NA0 = len(A0_LIST)
NWC = len(WC_LIST)
print(NA0,NWC,NPOL)
E = np.zeros(( 2, NA0, NWC)) # (x,m) x NPOL

sp.call(f"mkdir -p {DATA_DIR}", shell=True)

for A0_IND, A0 in enumerate( A0_LIST ):
    print(f"Reading {A0_IND+1} of {NA0}")
    for WC_IND, WC in enumerate( WC_LIST ):
        #print(f"\tReading {WC_IND+1} of {NWC}")
        A0 = round(A0,4)
        WC = round(WC,4)
        E[0,A0_IND,WC_IND] = np.loadtxt(f"oBr/QCHEM/PF_THETA_PHI/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat")[0] / 27.2114 * 630
        E[1,A0_IND,WC_IND] = np.loadtxt(f"mBr/QCHEM/PF_THETA_PHI/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat")[0] / 27.2114 * 630

######## PLOT 4 (2D) ########
# E_ORTHO - E_META

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, (E[0,:,:] - E[1,:,:]).T, kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1],500)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1],500)

X,Y = np.meshgrid(WC_FINE, A0_FINE)

divnorm=mpl.colors.TwoSlopeNorm(vcenter=0.0)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE).T, cmap="seismic", levels=1000, norm=divnorm )

np.savetxt(f"{DATA_DIR}/polariton_WCSCAN_A0SCAN_THETA_{THETA}_PHI_{PHI}.dat", f_xy(A0_FINE,WC_FINE))
np.savetxt(f"{DATA_DIR}/polariton_WCSCAN_A0SCAN_THETA_{THETA}_PHI_{PHI}_A0GRID.dat", A0_FINE)
np.savetxt(f"{DATA_DIR}/polariton_WCSCAN_A0SCAN_THETA_{THETA}_PHI_{PHI}_WCGRID.dat", WC_FINE)

plt.colorbar(pad=0.01)
#plt.legend()
#plt.xlim(A0_LIST[0],A0_LIST[-1])
#plt.ylim(0)
plt.xlabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.ylabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.savefig(f"{DATA_DIR}/polariton_WCSCAN_A0SCAN_THETA_{THETA}_PHI_{PHI}.jpg",dpi=600)
plt.clf()









V  = 0.15 * 10**3 / 0.529**3
WC = 1.8
A0 = np.sqrt( 2*np.pi / (WC/27.2114) / V )
print( "V, A0, E, = ",  0.15, A0, 514 * (WC/27.2114) * A0, "dE = ", f_xy(A0,WC) )
























# Cavity Volume vs. Cavity Frequency

NPTS = 500

f_xy = interp2d(A0_LIST, WC_LIST, (E[0,:,:] - E[1,:,:]).T, kind='cubic')
WC_FINE = np.linspace(1, 4, NPTS)

VOLUME_LIST  = np.linspace( 0.05,1.0,NPTS )  # nm**3
WV = np.outer( WC_FINE/27.2114, VOLUME_LIST * 10**3 / 0.529**3 )
A0 = np.sqrt( 2 * np.pi / WV ) # wc is first here !

f_xy_new = np.zeros( (NPTS,NPTS) )
for i in range(NPTS): # WC-LOOP
    for j in range(NPTS): # V-LOOP
        if ( A0[i,j] > 0.4 ): 
            f_xy_new[i,j] = float("Nan")
        else:
            f_xy_new[i,j] = f_xy(A0[i,j],WC_FINE[i])

W,V = np.meshgrid( WC_FINE, VOLUME_LIST )
divnorm=mpl.colors.TwoSlopeNorm(vcenter=0.0)
plt.contourf( W,V,f_xy_new.T, cmap="seismic",levels=1000, norm=divnorm )
plt.colorbar(pad=0.01)
plt.xlabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.ylabel("Cavity Volume, $\mathcal{V}$ (nm$^3$)",fontsize=15)
plt.savefig(f"{DATA_DIR}/polariton_VSCAN_WCCSCAN_THETA_{THETA}_PHI_{PHI}.jpg",dpi=600)
plt.clf()


np.savetxt(f"{DATA_DIR}/polariton_VSCAN_WCCSCAN_THETA_{THETA}_PHI_{PHI}_WGRID.dat", WC_FINE)
np.savetxt(f"{DATA_DIR}/polariton_VSCAN_WCCSCAN_THETA_{THETA}_PHI_{PHI}_VGRID.dat", VOLUME_LIST)
np.savetxt(f"{DATA_DIR}/polariton_VSCAN_WCCSCAN_THETA_{THETA}_PHI_{PHI}.dat",f_xy_new.T)






# Cavity Field Strength vs. Cavity Frequency

NPTS = 500

f_xy = interp2d(A0_LIST, WC_LIST, (E[0,:,:] - E[1,:,:]).T, kind='cubic')
WC_FINE = np.linspace(1, 4, NPTS)
f_xy_new = np.zeros( (NPTS,NPTS) )

FIELD_STRENGTH = np.linspace( 1,20,NPTS )  # V/nm
A0 = np.outer( 1/(WC_FINE/27.2114), FIELD_STRENGTH / 514 ) # A0 = E/WC

for i in range(NPTS): # WC-LOOP
    for j in range(NPTS): # F-LOOP
        if ( A0[i,j] > 0.4 ): 
            f_xy_new[i,j] = float("Nan")
        else:
            f_xy_new[i,j] = f_xy(A0[i,j],WC_FINE[i])

W,F = np.meshgrid( WC_FINE,  FIELD_STRENGTH )
divnorm=mpl.colors.TwoSlopeNorm(vcenter=0.0)
plt.contourf( W,F,f_xy_new.T, cmap="seismic",levels=1000, norm=divnorm )
plt.colorbar(pad=0.01)
plt.xlabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.ylabel("Cavity Field Strength, $|\mathcal{E}|$ (V/nm)",fontsize=15)
plt.savefig(f"{DATA_DIR}/polariton_FSCAN_WCCSCAN_THETA_{THETA}_PHI_{PHI}.jpg",dpi=600)
plt.clf()


np.savetxt(f"{DATA_DIR}/polariton_FSCAN_WCCSCAN_THETA_{THETA}_PHI_{PHI}_WGRID.dat", WC_FINE)
np.savetxt(f"{DATA_DIR}/polariton_FSCAN_WCCSCAN_THETA_{THETA}_PHI_{PHI}_FGRID.dat", FIELD_STRENGTH)
np.savetxt(f"{DATA_DIR}/polariton_FSCAN_WCCSCAN_THETA_{THETA}_PHI_{PHI}.dat",f_xy_new)








