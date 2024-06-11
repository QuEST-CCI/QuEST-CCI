import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import os
from scipy.interpolate import interp2d
import subprocess as sp

NF  = 5
NM  = 50
A0_LIST = np.arange( 0.0, 0.4+0.05, 0.05 )
#A0_LIST = np.array([ 0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
WC_LIST = np.arange( 0.0, 20.0+1.0, 1.0 )
#WC_LIST = np.array([ 0.0, 5.0, 10.0, 15.0, 20.0])
DATA_DIR = "PLOTS_DATA_ABS_ENERGY/"

NPOL = NM*NF
NA0 = len(A0_LIST)
NWC = len(WC_LIST)
print(NA0,NWC,NPOL)
E = np.zeros(( 4, NA0, NWC, NPOL )) # (o1,m1,m2,p2) x NPOL

sp.call(f"mkdir -p {DATA_DIR}", shell=True)

for A0_IND, A0 in enumerate( A0_LIST ):
    for WC_IND, WC in enumerate( WC_LIST ):
        A0 = round(A0,4)
        WC = round(WC,4)
        E[0,A0_IND,WC_IND,:] = np.loadtxt(f"oBr/QCHEM/PF_NF{NF}_NM{NM}_CORRECTED_DIPOLE/data_PF/E_1_0_0_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
        E[1,A0_IND,WC_IND,:] = np.loadtxt(f"mBr/QCHEM/PF_NF{NF}_NM{NM}_CORRECTED_DIPOLE/data_PF/E_1_0_0_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
        E[2,A0_IND,WC_IND,:] = np.loadtxt(f"mBr/QCHEM/PF_NF{NF}_NM{NM}_CORRECTED_DIPOLE/data_PF/E_0_1_0_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
        E[3,A0_IND,WC_IND,:] = np.loadtxt(f"pBr/QCHEM/PF_NF{NF}_NM{NM}_CORRECTED_DIPOLE/data_PF/E_0_1_0_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630

WCi = 2
EZERO = E[1,0,WCi,0]
plt.plot( A0_LIST, E[1,:,WCi,0] - EZERO, "-", c='black', label="Meta (X)" )
plt.plot( A0_LIST, E[0,:,WCi,0] - EZERO, "-", c='red', label="Ortho (X)" )
plt.plot( A0_LIST, E[2,:,WCi,0] - EZERO, "--",c='black',  label="Meta (Y)" )
plt.plot( A0_LIST, E[3,:,WCi,0] - EZERO, "--",c='red',  label="Para (Y)" )
plt.legend()
plt.title("$\omega_\mathrm{c}$ = %1.1f eV"%WC_LIST[WCi],fontsize=15)
plt.xlabel("Coupling Strength, $A_0$, (a.u.)",fontsize=15)
plt.ylabel("Energy (kcal/mol)",fontsize=15)
plt.savefig(f"{DATA_DIR}/Absolute_Energy.jpg", dpi=300)