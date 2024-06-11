import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import os
from scipy.interpolate import interp2d
import subprocess as sp

NF  = 5
NM  = 50
A0_LIST = np.array([0.3])
WC_LIST = np.array([2.5,10.0])
DATA_DIR = "PLOTS_DATA_XY_SCAN/"
sp.call(f"mkdir -p {DATA_DIR}", shell=True)

dA = 0.1
ANGLE_LIST = np.arange( 0, 2*np.pi+dA, dA )

NPOL   = NM*NF
NANGLE = len(ANGLE_LIST)
NA0    = len(A0_LIST)
NWC    = len(WC_LIST)
E = np.zeros(( 3, NANGLE, NA0, NWC, NPOL )) # (o,m,p) x NPOL

for ANGLEi,ANGLE in enumerate( ANGLE_LIST ): 
    for A0_IND, A0 in enumerate( A0_LIST ):
        #print(f"Reading {A0_IND+1} of {NA0}")
        for WC_IND, WC in enumerate( WC_LIST ):
            #print(f"\tReading {WC_IND+1} of {NWC}")
            A0    = round(A0,4)
            WC    = round(WC,4)
            ANGLE = round(ANGLE,2)
            E[0,ANGLEi,A0_IND,WC_IND,:] = np.loadtxt(f"oBr/QCHEM/PF_XY_SCAN/data_PF/E_{ANGLE}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
            E[1,ANGLEi,A0_IND,WC_IND,:] = np.loadtxt(f"mBr/QCHEM/PF_XY_SCAN/data_PF/E_{ANGLE}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
            E[2,ANGLEi,A0_IND,WC_IND,:] = np.loadtxt(f"pBr/QCHEM/PF_XY_SCAN/data_PF/E_{ANGLE}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630


# ORTHO
for WC_IND, WC in enumerate( WC_LIST ):
    for A0_IND, A0 in enumerate( A0_LIST ):
        plt.plot( ANGLE_LIST[:] * 180 / np.pi, E[0,:,A0_IND,WC_IND,0] - E[1,0,A0_IND,WC_IND,0], label="$A_0$"+f" = {round(A0_LIST[A0_IND],3)}" )
    plt.xlabel("Angle, $\\theta_{xy}$ (deg.)",fontsize=15)
    plt.ylabel("$E_0(Ortho)$ (kcal/mol)",fontsize=15)
    plt.title("$\omega_\mathrm{c}$ = "+f"{WC_LIST[WC_IND]} eV",fontsize=15)
    plt.xlim(0,180)
    plt.savefig(f"{DATA_DIR}/XY_SCAN_Ortho_WC_{WC_LIST[WC_IND]}.jpg", dpi=300)
    plt.clf()

# META
for WC_IND, WC in enumerate( WC_LIST ):
    for A0_IND, A0 in enumerate( A0_LIST ):
        plt.plot( ANGLE_LIST[:] * 180 / np.pi, E[1,:,A0_IND,WC_IND,0] - E[1,0,A0_IND,WC_IND,0], label="$A_0$"+f" = {round(A0_LIST[A0_IND],3)}" )
    plt.xlabel("Angle, $\\theta_{xy}$ (deg.)",fontsize=15)
    plt.ylabel("$E_0(Meta)$ (kcal/mol)",fontsize=15)
    plt.title("$\omega_\mathrm{c}$ = "+f"{WC_LIST[WC_IND]} eV",fontsize=15)
    plt.xlim(0,180)
    plt.savefig(f"{DATA_DIR}/XY_SCAN_Meta_WC_{WC_LIST[WC_IND]}.jpg", dpi=300)
    plt.clf()

# PARA
for WC_IND, WC in enumerate( WC_LIST ):
    for A0_IND, A0 in enumerate( A0_LIST ):
        plt.plot( ANGLE_LIST[:] * 180 / np.pi, E[1,:,A0_IND,WC_IND,0] - E[1,0,A0_IND,WC_IND,0], label="$A_0$"+f" = {round(A0_LIST[A0_IND],3)}" )
    plt.xlabel("Angle, $\\theta_{xy}$ (deg.)",fontsize=15)
    plt.ylabel("$E_0(Para)$ (kcal/mol)",fontsize=15)
    plt.title("$\omega_\mathrm{c}$ = "+f"{WC_LIST[WC_IND]} eV",fontsize=15)
    plt.xlim(0,180)
    plt.savefig(f"{DATA_DIR}/XY_SCAN_Para_WC_{WC_LIST[WC_IND]}.jpg", dpi=300)
    plt.clf()




# ORTHO/META
for WC_IND, WC in enumerate( WC_LIST ):
    for A0_IND, A0 in enumerate( A0_LIST ):
        plt.plot( ANGLE_LIST[:] * 180 / np.pi, E[0,:,A0_IND,WC_IND,0] - E[1,:,A0_IND,WC_IND,0], label="$A_0$"+f" = {round(A0_LIST[A0_IND],3)}" )
    plt.xlabel("Angle, $\\theta_{xy}$ (deg.)",fontsize=15)
    plt.ylabel("$E_0(Ortho) - E_0(Meta)$ (kcal/mol)",fontsize=15)
    plt.title("$\omega_\mathrm{c}$ = "+f"{WC_LIST[WC_IND]} eV",fontsize=15)
    plt.xlim(0,180)
    plt.savefig(f"{DATA_DIR}/XY_SCAN_OrthoMeta_WC_{WC_LIST[WC_IND]}.jpg", dpi=300)
    plt.clf()


# PARA/META
for WC_IND, WC in enumerate( WC_LIST ):
    for A0_IND, A0 in enumerate( A0_LIST ):
        plt.plot( ANGLE_LIST[:] * 180 / np.pi, E[2,:,A0_IND,WC_IND,0] - E[1,:,A0_IND,WC_IND,0], label="$A_0$"+f" = {round(A0_LIST[A0_IND],3)}" )
    plt.xlabel("Angle, $\\theta_{xy}$ (deg.)",fontsize=15)
    plt.ylabel("$E_0(Para) - E_0(Meta)$ (kcal/mol)",fontsize=15)
    plt.title("$\omega_\mathrm{c}$ = "+f"{WC_LIST[WC_IND]} eV",fontsize=15)
    plt.xlim(0,180)
    plt.savefig(f"{DATA_DIR}/XY_SCAN_ParaMeta_WC_{WC_LIST[WC_IND]}.jpg", dpi=300)
    plt.clf()