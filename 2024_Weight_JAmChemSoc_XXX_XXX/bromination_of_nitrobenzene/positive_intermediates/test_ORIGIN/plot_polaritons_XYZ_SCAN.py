import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import os
from scipy.interpolate import interp2d
import subprocess as sp

COHERENT_STATE = ["FOCK_BASIS", "COHERENT_BASIS"]
ORIGIN = ["ORIGIN", "SHIFTED"]

NF  = 5
NM  = 15
A0_LIST = np.array([0.1,0.2,0.3])
WC_LIST = np.array([1.8])

dA = 0.2
THETA_LIST = np.arange( 0, np.pi+dA, dA )
PHI_LIST   = np.arange( 0, 2*np.pi+dA, dA )

NPOL   = NM*NF
NTHETA = len(THETA_LIST)
NPHI   = len(PHI_LIST)
NA0    = len(A0_LIST)
NWC    = len(WC_LIST)
E = np.zeros(( 2, 2, NTHETA, NPHI, NA0, NWC, NPOL ))

DATA_DIR = "PLOTS/"
sp.call(f"mkdir -p {DATA_DIR}", shell=True)

for origi,orig in enumerate( ORIGIN ):
    for BASISi,BASIS in enumerate( COHERENT_STATE ):
        for THETAi,THETA in enumerate( THETA_LIST ):
            print(f"Reading {THETAi+1} of {NTHETA}")
            for PHIi,PHI in enumerate( PHI_LIST ): 
                for A0_IND, A0 in enumerate( A0_LIST ):
                    #print(f"Reading {A0_IND+1} of {NA0}")
                    for WC_IND, WC in enumerate( WC_LIST ):
                        #print(f"\tReading {WC_IND+1} of {NWC}")
                        A0    = round(A0,4)
                        WC    = round(WC,4)
                        THETA = round(THETA,2)
                        PHI   = round(PHI,2)
                        E[BASISi,origi,THETAi,PHIi,A0_IND,WC_IND,:] = np.loadtxt(f"{BASIS}/{orig}/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
                        

EZERO = E[1,1,len(THETA_LIST)//2,len(PHI_LIST)//2,0,len(WC_LIST)//2,0]
plt.plot( THETA_LIST, E[0,0,:,len(PHI_LIST)//2,1,len(WC_LIST)//2,0] - EZERO, "-",  label=f"FOCK, ORIGIN" )
#plt.plot( THETA_LIST, E[0,1,:,len(PHI_LIST)//2,1,len(WC_LIST)//2,0] - EZERO, "--",  label=f"FOCK, SHIFTED" )
plt.plot( THETA_LIST, E[1,0,:,len(PHI_LIST)//2,1,len(WC_LIST)//2,0] - EZERO, "-o", label=f"COHERENT, ORIGIN" )
plt.plot( THETA_LIST, E[1,1,:,len(PHI_LIST)//2,1,len(WC_LIST)//2,0] - EZERO, "-.", label=f"COHERENT, SHIFTED" )
plt.legend()
plt.savefig(f"{DATA_DIR}/COMPARISON.jpg", dpi=300)
plt.clf()

plt.plot( THETA_LIST, E[0,0,:,len(PHI_LIST)//2,1,len(WC_LIST)//2,0] - E[1,1,:,len(PHI_LIST)//2,1,len(WC_LIST)//2,0], "-",  label=f"FOCK, ORIGIN" )
#plt.plot( THETA_LIST, E[0,1,:,len(PHI_LIST)//2,1,len(WC_LIST)//2,0] - E[1,1,:,len(PHI_LIST)//2,1,len(WC_LIST)//2,0], "--",  label=f"FOCK, SHIFTED" )
plt.plot( THETA_LIST, E[1,0,:,len(PHI_LIST)//2,1,len(WC_LIST)//2,0] - E[1,1,:,len(PHI_LIST)//2,1,len(WC_LIST)//2,0], "-o", label=f"COHERENT, ORIGIN" )
plt.plot( THETA_LIST, E[1,1,:,len(PHI_LIST)//2,1,len(WC_LIST)//2,0] - E[1,1,:,len(PHI_LIST)//2,1,len(WC_LIST)//2,0], "-.", label=f"COHERENT, SHIFTED" )
plt.legend()
plt.savefig(f"{DATA_DIR}/DIFFERENCE.jpg", dpi=300)
plt.clf()

