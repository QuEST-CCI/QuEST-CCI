import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import os
from scipy.interpolate import interp2d
import subprocess as sp

NF  = 5
NM  = 50
#A0_LIST = np.array([0, 0.4+0.05, 0.05])
A0_LIST = np.arange(0.0, 0.2+0.05, 0.05)
#WC_LIST = np.array([0.0, 5.0, 10.0, 15.0, 20.0])
WC_LIST = np.arange(0.0,20.0+0.1,0.1)
DATA_DIR = "PLOTS_DATA_XYZ_SCAN_AVERAGE/"
sp.call(f"mkdir -p {DATA_DIR}", shell=True)

dA = 0.4
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
                E[0,THETAi,PHIi,A0_IND,WC_IND,:] = np.loadtxt(f"oBr/QCHEM/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
                E[1,THETAi,PHIi,A0_IND,WC_IND,:] = np.loadtxt(f"mBr/QCHEM/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
                E[2,THETAi,PHIi,A0_IND,WC_IND,:] = np.loadtxt(f"pBr/QCHEM/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630

import matplotlib as mpl
cmap=mpl.colors.LinearSegmentedColormap.from_list('rg',[ "darkblue", "blue", "white", "darkred", "red" ], N=256)

"""
# PARA/META
for WC_IND, WC in enumerate( WC_LIST ):
    for A0_IND, A0 in enumerate( A0_LIST ):
        A0    = round(A0,4)
        WC    = round(WC,4)
        DATA = E[2,:,:,A0_IND,WC_IND,0] - E[1,:,:,A0_IND,WC_IND,0]
        VMIN = np.min( DATA )
        VMAX = np.max( DATA )
        divnorm=mpl.colors.TwoSlopeNorm(vcenter=0.0)
        THETA,PHI = np.meshgrid( THETA_LIST* 180 / np.pi, PHI_LIST* 180 / np.pi )
        plt.contourf( PHI, THETA, DATA.T, cmap=cmap, levels=1000, norm=divnorm )
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
            
# ORTHO/META
for WC_IND, WC in enumerate( WC_LIST ):
    for A0_IND, A0 in enumerate( A0_LIST ):
        A0    = round(A0,4)
        WC    = round(WC,4)
        DATA = E[0,:,:,A0_IND,WC_IND,0] - E[1,:,:,A0_IND,WC_IND,0]
        VMIN = np.min( DATA )
        VMAX = np.max( DATA )
        divnorm=mpl.colors.TwoSlopeNorm(vcenter=0.0)
        THETA,PHI = np.meshgrid( THETA_LIST* 180 / np.pi, PHI_LIST* 180 / np.pi )
        plt.contourf( PHI, THETA, DATA.T, cmap=cmap, levels=1000, norm=divnorm )
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
"""


### COMPUTE AVERAGE OVER ORIENTATIONS ###
AVERAGE = np.zeros(( 2, NA0, NWC )) # o-m, p-m
AVERAGE[0,:,:] = np.average( E[0,:,:,:,:,0] - E[1,:,:,:,:,0], axis=(0,1) )
AVERAGE[1,:,:] = np.average( E[2,:,:,:,:,0] - E[1,:,:,:,:,0], axis=(0,1) )

print(AVERAGE[0,:,:])
print(AVERAGE[1,:,:])

print( AVERAGE.shape )
print( A0_LIST.shape )
print( WC_LIST.shape )

A0_FINE = np.linspace( A0_LIST[0], A0_LIST[-1], 500 )
WC_FINE = np.linspace( WC_LIST[0], WC_LIST[-1], 500 )
#f_om    = interp2d( WC_LIST, A0_LIST, AVERAGE[0,:,:], kind='cubic' )
#f_pm    = interp2d( WC_LIST, A0_LIST, AVERAGE[1,:,:], kind='cubic' )

A,W = np.meshgrid( A0_LIST, WC_LIST )
VMIN = np.min( AVERAGE[0,:,:] )
VMAX = np.max( AVERAGE[0,:,:] )
divnorm=mpl.colors.TwoSlopeNorm(vcenter=0.0)
#plt.contourf( A,W,f_om(A0_FINE,WC_FINE), cmap=cmap, norm=divnorm )
plt.imshow( AVERAGE[0,:,:], cmap=cmap, norm=divnorm )
plt.xlabel("Coupling Strength, $A_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_\mathrm{c}$ (eV)",fontsize=15)
#plt.xlim(0,A0_LIST[-1])
#plt.ylim(0,WC_LIST[-1])
cbar = plt.colorbar(pad=0.01)
cbar.set_label("$E_0(Ortho) - E_0(Meta)$ (kcal/mol)", size=15)
cbar.ax.tick_params(labelsize=10 ) 
plt.savefig(f"{DATA_DIR}/XYZ_SCAN_OrthoMeta_AVERAGE.jpg", dpi=300)
plt.clf()

A,W = np.meshgrid( A0_LIST, WC_LIST )
VMIN = np.min( AVERAGE[0,:,:] )
VMAX = np.max( AVERAGE[0,:,:] )
divnorm=mpl.colors.TwoSlopeNorm(vcenter=0.0)
#plt.contourf( A,W,f_pm(A0_FINE,WC_FINE), cmap=cmap, norm=divnorm )
plt.imshow( AVERAGE[1,:,:], cmap=cmap, norm=divnorm )
plt.xlabel("Coupling Strength, $A_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_\mathrm{c}$ (eV)",fontsize=15)
#plt.xlim(0,A0_LIST[-1])
#plt.ylim(0,WC_LIST[-1])
cbar = plt.colorbar(pad=0.01)
cbar.set_label("$E_0(Para) - E_0(Meta)$ (kcal/mol)", size=15)
cbar.ax.tick_params(labelsize=10 ) 
plt.savefig(f"{DATA_DIR}/XYZ_SCAN_ParaMeta_AVERAGE.jpg", dpi=300)
plt.clf()