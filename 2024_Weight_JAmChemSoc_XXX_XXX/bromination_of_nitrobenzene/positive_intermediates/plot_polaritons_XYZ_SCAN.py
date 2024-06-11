import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import os
from scipy.interpolate import interp2d
import subprocess as sp

COHERENT_STATE = True

NF  = 5
NM  = 50
A0_LIST = np.array([0.3])
WC_LIST = np.array([1.8])
DATA_DIR = f"PLOTS_DATA_XYZ_SCAN_COHERENT_STATE_{COHERENT_STATE}/"
sp.call(f"mkdir -p {DATA_DIR}", shell=True)

dA = 0.2
THETA_LIST = np.arange( 0, np.pi+dA, dA )
PHI_LIST   = np.arange( 0, 2*np.pi+dA, dA )

NPOL   = NM*NF
NTHETA = len(THETA_LIST)
NPHI   = len(PHI_LIST)
NA0    = len(A0_LIST)
NWC    = len(WC_LIST)
E = np.zeros(( 3, NTHETA, NPHI, NA0, NWC, NPOL )) # (o,m,p) x NPOL

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
                if ( COHERENT_STATE == False ):
                    E[0,THETAi,PHIi,A0_IND,WC_IND,:] = np.loadtxt(f"oBr/QCHEM/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
                    E[1,THETAi,PHIi,A0_IND,WC_IND,:] = np.loadtxt(f"mBr/QCHEM/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
                    E[2,THETAi,PHIi,A0_IND,WC_IND,:] = np.loadtxt(f"pBr/QCHEM/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
                else:
                    E[0,THETAi,PHIi,A0_IND,WC_IND,:] = np.loadtxt(f"oBr/QCHEM/PF_COHERENT_STATE/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
                    E[1,THETAi,PHIi,A0_IND,WC_IND,:] = np.loadtxt(f"mBr/QCHEM/PF_COHERENT_STATE/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
                    E[2,THETAi,PHIi,A0_IND,WC_IND,:] = np.loadtxt(f"pBr/QCHEM/PF_COHERENT_STATE/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
                   

import matplotlib as mpl
cmap=mpl.colors.LinearSegmentedColormap.from_list('rg',[ "darkblue", "blue", "white", "darkred", "red" ], N=256)

THETA_FINE = np.linspace( THETA_LIST[0], THETA_LIST[-1], 500 )
PHI_FINE   = np.linspace( PHI_LIST[0], PHI_LIST[-1], 500 )

# PARA/META
for WC_IND, WC in enumerate( WC_LIST ):
    for A0_IND, A0 in enumerate( A0_LIST ):
        A0    = round(A0,4)
        WC    = round(WC,4)
        DATA = E[2,:,:,A0_IND,WC_IND,0] - E[1,:,:,A0_IND,WC_IND,0]
        f = interp2d( THETA_LIST, PHI_LIST, DATA.T, kind='cubic' )
        if ( np.max(DATA) < 0 or np.min(DATA) > 0 ):
            VCENTER = 0.5 * (np.min(DATA) + np.max(DATA))
        else:
            VCENTER = 0.0
        divnorm=mpl.colors.TwoSlopeNorm(vcenter=VCENTER)
        THETA,PHI = np.meshgrid( THETA_FINE* 180 / np.pi, PHI_FINE* 180 / np.pi )
        
        print( "Para - Meta", np.min( f(THETA_FINE,PHI_FINE) ) )
        
        plt.contourf( PHI, THETA, f(THETA_FINE,PHI_FINE), cmap=cmap, levels=1000, norm=divnorm )
        plt.xlabel("Angle, $\phi$ (deg.)",fontsize=15)
        plt.ylabel("Angle, $\\theta$ (deg.)",fontsize=15)
        plt.title("$A_0$"+f" = {A0} a.u.; " + "$\omega_\mathrm{c}$ = "+f"{WC} eV",fontsize=15)
        plt.xlim(0,360)
        plt.ylim(0,180)
        cbar = plt.colorbar(pad=0.01)
        cbar.set_label("$E_0(Para) - E_0(Meta)$ (kcal/mol)", size=15)
        cbar.ax.tick_params(labelsize=10 ) 
        plt.savefig(f"{DATA_DIR}/XYZ_SCAN_ParaMeta_A0_{A0}_WC_{WC}.jpg", dpi=300)
        plt.clf()
            
        np.savetxt(f"{DATA_DIR}/XYZ_SCAN_ParaMeta_A0_{A0}_WC_{WC}.dat", f(THETA_FINE,PHI_FINE), fmt="%1.4f" )
        np.savetxt(f"{DATA_DIR}/XYZ_SCAN_ParaMeta_A0_{A0}_WC_{WC}_THETA.dat", THETA_FINE, fmt="%1.4f" )
        np.savetxt(f"{DATA_DIR}/XYZ_SCAN_ParaMeta_A0_{A0}_WC_{WC}_PHI.dat", PHI_FINE, fmt="%1.4f" )


# ORTHO/META
for WC_IND, WC in enumerate( WC_LIST ):
    for A0_IND, A0 in enumerate( A0_LIST ):
        A0    = round(A0,4)
        WC    = round(WC,4)
        DATA = E[0,:,:,A0_IND,WC_IND,0] - E[1,:,:,A0_IND,WC_IND,0]
        f = interp2d( THETA_LIST, PHI_LIST, DATA.T, kind='cubic' )
        if ( np.max(DATA) < 0 or np.min(DATA) > 0 ):
            VCENTER = 0.5 * (np.min(DATA) + np.max(DATA))
        else:
            VCENTER = 0.0
        divnorm=mpl.colors.TwoSlopeNorm(vcenter=VCENTER)
        #THETA,PHI = np.meshgrid( THETA_LIST* 180 / np.pi, PHI_LIST* 180 / np.pi )
        THETA,PHI = np.meshgrid( THETA_FINE* 180 / np.pi, PHI_FINE* 180 / np.pi )
        
        print( "Ortho - Meta", np.min( f(THETA_FINE,PHI_FINE) ) )
        
        plt.contourf( PHI, THETA, f(THETA_FINE,PHI_FINE), cmap=cmap, levels=1000, norm=divnorm )
        plt.xlabel("Angle, $\phi$ (deg.)",fontsize=15)
        plt.ylabel("Angle, $\\theta$ (deg.)",fontsize=15)
        plt.title("$A_0$"+f" = {A0} a.u.; " + "$\omega_\mathrm{c}$ = "+f"{WC} eV",fontsize=15)
        plt.xlim(0,360)
        plt.ylim(0,180)
        cbar = plt.colorbar(pad=0.01)
        cbar.set_label("$E_0(Ortho) - E_0(Meta)$ (kcal/mol)", size=15)
        cbar.ax.tick_params(labelsize=10 ) 
        plt.savefig(f"{DATA_DIR}/XYZ_SCAN_OrthoMeta_A0_{A0}_WC_{WC}.jpg", dpi=300)
        plt.clf()

        np.savetxt(f"{DATA_DIR}/XYZ_SCAN_OrthoMeta_A0_{A0}_WC_{WC}.dat", f(THETA_FINE,PHI_FINE), fmt="%1.4f" )
        np.savetxt(f"{DATA_DIR}/XYZ_SCAN_OrthoMeta_A0_{A0}_WC_{WC}_THETA.dat", THETA_FINE * 180 / np.pi, fmt="%1.4f" )
        np.savetxt(f"{DATA_DIR}/XYZ_SCAN_OrthoMeta_A0_{A0}_WC_{WC}_PHI.dat", PHI_FINE * 180 / np.pi, fmt="%1.4f" )

