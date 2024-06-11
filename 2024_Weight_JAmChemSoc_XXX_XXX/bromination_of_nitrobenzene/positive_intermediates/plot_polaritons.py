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
WC_LIST = np.arange( 0.0, 5.0+0.1, 0.1 )
#WC_LIST = np.array([ 0.0, 5.0, 10.0, 15.0, 20.0])
DATA_DIR = "PLOTS_DATA/"

EVEC_INTS = np.array([ 0,1,0 ])
EVEC_NORM = EVEC_INTS / np.linalg.norm(EVEC_INTS)
EVEC_OUT = "_".join(map(str,EVEC_INTS))


NPOL = NM*NF
NA0 = len(A0_LIST)
NWC = len(WC_LIST)
print(NA0,NWC,NPOL)
E = np.zeros(( 3, NA0, NWC, NPOL )) # (o,m,p) x NPOL
PHOT_0 = np.zeros(( 3, NA0, NWC )) # (o,m,p) x NPOL
EXCI_0 = np.zeros(( 3, NA0, NWC )) # (o,m,p) x NPOL
#U = np.zeros(( 2, NA0, NWC, NPOL, NPOL )) # (o,m,p) x NPOL
MU       = np.zeros(( 3, NM, NM ))
MU_POL_0 = np.zeros(( 3, NA0, NWC ))

sp.call(f"mkdir -p {DATA_DIR}", shell=True)

a_op = np.zeros((NF,NF))
for m in range(1,NF):
    a_op[m,m-1] = np.sqrt(m)
a_op = a_op.T
N_a_OP = a_op.T @ a_op

c_op = np.zeros((NM,NM))
for m in range(1,NM):
    c_op[m,m-1] = np.sqrt(m)
c_op = c_op.T
N_c_OP = c_op.T @ c_op

I_m = np.identity(NM)
I_ph = np.identity(NF)
PHOT_OP_FULL = np.kron( I_m, N_a_OP )
EXCI_OP_FULL = np.kron( N_c_OP, I_ph  )

#if ( not os.path.isfile(f"{DATA_DIR}/E_POL.npy") ):
if ( True ):

    TMP = np.load(f"oBr/QCHEM/PLOTS_DATA/DIPOLE_RPA.dat.npy")[:NM,:NM]
    MU[0,:,:] = np.einsum( "JKd,d->JK", TMP[:,:,:], np.array([1,0,0]) )

    TMP = np.load(f"mBr/QCHEM/PLOTS_DATA/DIPOLE_RPA.dat.npy")[:NM,:NM]
    MU[1,:,:] = np.einsum( "JKd,d->JK", TMP[:,:,:], np.array([1,0,0]) )
    
    TMP = np.load(f"pBr/QCHEM/PLOTS_DATA/DIPOLE_RPA.dat.npy")[:NM,:NM]
    MU[2,:,:] = np.einsum( "JKd,d->JK", TMP[:,:,:], np.array([1,0,0]) )
    
    np.save(f"{DATA_DIR}/MU.npy", MU)

    MU_POL_OP = np.array([np.kron( MU[0,:,:], I_ph ), \
                          np.kron( MU[1,:,:], I_ph ), \
                          np.kron( MU[2,:,:], I_ph ) ])

    for A0_IND, A0 in enumerate( A0_LIST ):
        #print(f"Reading {A0_IND+1} of {NA0}")
        for WC_IND, WC in enumerate( WC_LIST ):
            #print(f"\tReading {WC_IND+1} of {NWC}")
            A0 = round(A0,4)
            WC = round(WC,4)
            E[0,A0_IND,WC_IND,:] = np.loadtxt(f"oBr/QCHEM/PF_NF{NF}_NM{NM}_CORRECTED_DIPOLE/data_PF/E_{EVEC_OUT}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
            E[1,A0_IND,WC_IND,:] = np.loadtxt(f"mBr/QCHEM/PF_NF{NF}_NM{NM}_CORRECTED_DIPOLE/data_PF/E_{EVEC_OUT}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
            E[2,A0_IND,WC_IND,:] = np.loadtxt(f"pBr/QCHEM/PF_NF{NF}_NM{NM}_CORRECTED_DIPOLE/data_PF/E_{EVEC_OUT}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
            TMP = np.load(f"oBr/QCHEM/PF_NF{NF}_NM{NM}_CORRECTED_DIPOLE/data_PF/U_{EVEC_OUT}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat.npy")[:,0]
            PHOT_0[0,A0_IND,WC_IND]   = TMP @ PHOT_OP_FULL[:,:] @ TMP
            EXCI_0[0,A0_IND,WC_IND]   = TMP @ EXCI_OP_FULL[:,:] @ TMP
            MU_POL_0[0,A0_IND,WC_IND] = TMP @ MU_POL_OP[0,:,:]  @ TMP
            TMP = np.load(f"mBr/QCHEM/PF_NF{NF}_NM{NM}_CORRECTED_DIPOLE/data_PF/U_{EVEC_OUT}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat.npy")[:,0]
            PHOT_0[1,A0_IND,WC_IND]   = TMP @ PHOT_OP_FULL[:,:] @ TMP
            EXCI_0[1,A0_IND,WC_IND]   = TMP @ EXCI_OP_FULL[:,:] @ TMP
            MU_POL_0[1,A0_IND,WC_IND] = TMP @ MU_POL_OP[1,:,:]  @ TMP
            TMP = np.load(f"pBr/QCHEM/PF_NF{NF}_NM{NM}_CORRECTED_DIPOLE/data_PF/U_{EVEC_OUT}_A0_{round(A0,6)}_WC_{round(WC,6)}_NF_{NF}_NM_{NM}.dat.npy")[:,0]
            PHOT_0[2,A0_IND,WC_IND]   = TMP @ PHOT_OP_FULL[:,:] @ TMP
            EXCI_0[2,A0_IND,WC_IND]   = TMP @ EXCI_OP_FULL[:,:] @ TMP
            MU_POL_0[2,A0_IND,WC_IND] = TMP @ MU_POL_OP[2,:,:]  @ TMP

            # print(f"\tParameters (A0,WC): ({A0},{WC})")
            # print("\t\tE(oBr):", \
            #                  np.round(E[0,A0_IND,WC_IND,0]-E[0,0,0,0],5), np.round(E[1,A0_IND,WC_IND,0]-E[1,0,0,0],5), \
            #                  np.round(E[0,A0_IND,WC_IND,0]-E[1,A0_IND,WC_IND,0],5) )
            # print("\t\tPHOT:", 
            #                  np.round(PHOT_0[0,A0_IND,WC_IND],5), np.round(EXCI_0[0,A0_IND,WC_IND],5), \
            #                  np.round(PHOT_0[1,A0_IND,WC_IND],5), np.round(EXCI_0[1,A0_IND,WC_IND],5), \
            #                  np.round(PHOT_0[2,A0_IND,WC_IND],5), np.round(EXCI_0[2,A0_IND,WC_IND],5) )
            # print( "\t\tDIP:", 
            #                  np.round(MU_POL_0[0,A0_IND,WC_IND],5), np.round(MU[0,0,0],5), \
            #                  np.round(MU_POL_0[1,A0_IND,WC_IND],5), np.round(MU[1,0,0],5), \
            #                  np.round(MU_POL_0[2,A0_IND,WC_IND],5), np.round(MU[2,0,0],5) )
    
else:
    print("Found data. Reading it.")
    E      = np.load(f"{DATA_DIR}/E_POL.npy")
    #U      = np.load(f"{DATA_DIR}/U_POL.npy")
    PHOT_0 = np.load(f"{DATA_DIR}/PHOT_0_POL.npy")
    EXCI_0 = np.load(f"{DATA_DIR}/EXCI_0_POL.npy")
    MU        = np.load(f"{DATA_DIR}/MU.npy")
    MU_POL_0  = np.load(f"{DATA_DIR}/MU_POL_0.npy")

"""


######## PLOT 5 (2D) DIP oBr ########

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, np.abs(MU_POL_0[0,:,:].T), kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1], 500)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1], 500)

X,Y = np.meshgrid(A0_FINE,WC_FINE)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="hot", levels=500)

np.savetxt(f"{DATA_DIR}/2D_CONTOUR_DIP_oBr_{EVEC_OUT}.dat", f_xy(A0_LIST,WC_LIST))

plt.colorbar(pad=0.01)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.title("oBr Ground State Dipole, $|\mu_{00}| (a.u.)$",fontsize=15)
plt.savefig(f"{DATA_DIR}/DIP_WCSCAN_A0SCAN_oBr_{EVEC_OUT}.jpg",dpi=600)
plt.clf()


######## PLOT 5 (2D) DIP mBr ########

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, np.abs(MU_POL_0[1,:,:].T), kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1], 500)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1], 500)

X,Y = np.meshgrid(A0_FINE,WC_FINE)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="hot", levels=500)

np.savetxt(f"{DATA_DIR}/2D_CONTOUR_DIP_mBr_{EVEC_OUT}.dat", f_xy(A0_LIST,WC_LIST))

plt.colorbar(pad=0.01)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.title("mBr Ground State Dipole, $|\mu_{00}|$ (a.u.)",fontsize=15)
plt.savefig(f"{DATA_DIR}/DIP_WCSCAN_A0SCAN_mBr_{EVEC_OUT}.jpg",dpi=600)
plt.clf()




######## PLOT 5 (2D) PHOT oBr ########

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, PHOT_0[0,:,:].T, kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1], 100)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1], 100)

X,Y = np.meshgrid(A0_FINE,WC_FINE)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="hot_r", levels=1000)

np.savetxt(f"{DATA_DIR}/2D_CONTOUR_PHOT_oBr_{EVEC_OUT}.dat", f_xy(A0_LIST,WC_LIST))

plt.colorbar(pad=0.01)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.title("oBr Average Photon Number, $\hat{a}^\dag \hat{a}$",fontsize=15)
plt.savefig(f"{DATA_DIR}/PHOT_WCSCAN_A0SCAN_oBr_{EVEC_OUT}.jpg",dpi=600)
plt.clf()

######## PLOT 5 (2D) PHOT mBr ########

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, PHOT_0[1,:,:].T, kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1], 100)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1], 100)

X,Y = np.meshgrid(A0_FINE,WC_FINE)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="hot_r", levels=1000)

np.savetxt(f"{DATA_DIR}/2D_CONTOUR_PHOT_mBr_{EVEC_OUT}.dat", f_xy(A0_LIST,WC_LIST) )

plt.colorbar(pad=0.01)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.title("mBr Average Photon Number, $\hat{a}^\dag \hat{a}$",fontsize=15)
plt.savefig(f"{DATA_DIR}/PHOT_WCSCAN_A0SCAN_mBr_{EVEC_OUT}.jpg",dpi=600)
plt.clf()

######## PLOT 5 (2D) EXCI oBr ########

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, EXCI_0[0,:,:].T, kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1], 100)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1], 100)

X,Y = np.meshgrid(A0_FINE,WC_FINE)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="hot_r", levels=1000)

np.savetxt(f"{DATA_DIR}/2D_CONTOUR_EXCI_oBr_{EVEC_OUT}.dat", f_xy(A0_LIST,WC_LIST))

plt.colorbar(pad=0.01)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.title("oBr Average Exciton Number, $\hat{c}^\dag \hat{c}$",fontsize=15)
plt.savefig(f"{DATA_DIR}/EXCI_WCSCAN_A0SCAN_oBr_{EVEC_OUT}.jpg",dpi=600)
plt.clf()

######## PLOT 5 (2D) EXCI mBr ########

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, EXCI_0[1,:,:].T, kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1], 100)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1], 100)

X,Y = np.meshgrid(A0_FINE,WC_FINE)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="hot_r", levels=1000)

np.savetxt(f"{DATA_DIR}/2D_CONTOUR_EXCI_mBr_{EVEC_OUT}.dat", f_xy(A0_LIST,WC_LIST))

plt.colorbar(pad=0.01)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.title("mBr Average Exciton Number, $\hat{c}^\dag \hat{c}$",fontsize=15)
plt.savefig(f"{DATA_DIR}/EXCI_WCSCAN_A0SCAN_mBr_{EVEC_OUT}.jpg",dpi=600)
plt.clf()
















######## PLOT 1 ########
WC_IND =len(WC_LIST)//2
EREF = E[1,0,WC_IND,0] # E_{META} ( A0 = 0 )
plt.plot( A0_LIST, E[0,:,WC_IND,0] - EREF, "-", lw=4, c="black", label="Ortho ($\Phi_0$)" )
plt.plot( A0_LIST, E[1,:,WC_IND,0] - EREF, "-", lw=4, c="red", label="Meta ($\Phi_0$)" )
plt.plot( A0_LIST, E[2,:,WC_IND,0] - EREF, "-", lw=4, c="blue", label="Para ($\Phi_0$)" )
#plt.plot( A0_LIST, E[0,:,WC_IND,1] - EREF, ".", lw=4, c="black", label="$\Phi_1$" )
plt.plot( A0_LIST, E[1,:,WC_IND,1] - EREF, ".", lw=4, c="red")#, label="Meta ($\Phi_1$)" )
plt.plot( A0_LIST, E[2,:,WC_IND,1] - EREF, ".", lw=4, c="blue")#, label="Para ($\Phi_1$)" )
#plt.plot( A0_LIST[::3], E[2,::3,WC_IND,2] - EREF, "o", lw=4, c="blue")#, label="$\Phi_2$" )

plt.legend()
plt.xlim(A0_LIST[0],A0_LIST[-1])
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("E$^{G.S.}$(A$_0$) - E$^{G.S.}_{Ortho}$(A$_0$ = 0.0) (kcal/mol)",fontsize=15)
plt.title("$\omega_c$ = %2.2f eV" % (WC_LIST[WC_IND]),fontsize=15)
plt.savefig(f"{DATA_DIR}/polariton_A0SCAN_WC{round(WC_LIST[WC_IND],4)}eV_{EVEC_OUT}.jpg",dpi=600)
plt.clf()


######## PLOT 2 ########
plt.plot( A0_LIST, E[0,:,WC_IND,0] - E[1,0,WC_IND,0], lw=4, label="Ortho" )
plt.plot( A0_LIST, E[1,:,WC_IND,0] - E[1,0,WC_IND,0], lw=4, label="Meta" )
plt.plot( A0_LIST, E[2,:,WC_IND,0] - E[1,0,WC_IND,0], lw=4, label="Para" )

plt.legend()
plt.xlim(A0_LIST[0],A0_LIST[-1])
#plt.ylim(0,15)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("E$^{G.S.}_j$(A$_0$) - E$^{G.S.}_j$(A$_0$ = 0) (kcal/mol)",fontsize=15)
plt.savefig(f"{DATA_DIR}/polariton_shift_A0SCAN_WC{round(WC_LIST[WC_IND],4)}eV_{EVEC_OUT}.jpg",dpi=600)
plt.clf()

######## PLOT 3 ########
plt.plot( WC_LIST, E[0,NA0//2,:,0] - E[1,NA0//2,0,0], lw=4, label="Ortho" )
plt.plot( WC_LIST, E[1,NA0//2,:,0] - E[1,NA0//2,0,0], lw=4, label="Meta" )
plt.plot( WC_LIST, E[2,NA0//2,:,0] - E[1,NA0//2,0,0], lw=4, label="Para" )

plt.legend()
plt.xlim(WC_LIST[0],WC_LIST[-1])
#plt.ylim(0)
plt.xlabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.ylabel("E$^{G.S.}_j$($\omega_c$) - E$^{G.S.}_{Ortho}$($\omega_c$ = 0.0 eV) (kcal/mol)",fontsize=15)
plt.savefig(f"{DATA_DIR}/polariton_shift_A0{round(A0_LIST[NA0//2],4)}_WCSCAN_{EVEC_OUT}.jpg",dpi=600)
plt.clf()
"""


######## PLOT 4 (2D) ########
# E_ORTHO - E_META

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, (E[0,:,:,0] - E[1,:,:,0]).T, kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1],100)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1],100)

X,Y = np.meshgrid(A0_FINE,WC_FINE)

VMAX = np.max( np.abs(f_xy(A0_FINE,WC_FINE)) )
VMIN = - VMAX

print("Ortho - Meta", f_xy(0.31,1.8) )

plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="seismic", levels=500 ) # , norm=Normalize(vmin=VMIN, vmax=VMAX)

np.savetxt(f"{DATA_DIR}/2D_CONTOUR_{EVEC_OUT}_01.dat", f_xy(A0_FINE,WC_FINE))

plt.colorbar(pad=0.01)
#plt.legend()
#plt.xlim(A0_LIST[0],A0_LIST[-1])
#plt.ylim(0)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.savefig(f"{DATA_DIR}/polariton_WCSCAN_A0SCAN_{EVEC_OUT}_01.jpg",dpi=600)
plt.clf()





######## PLOT 4 (2D) ########
# E_ORTHO - E_META

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, (E[2,:,:,0] - E[1,:,:,0]).T, kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1],100)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1],100)

X,Y = np.meshgrid(A0_FINE,WC_FINE)

VMAX = np.max( np.abs(f_xy(A0_FINE,WC_FINE)) )
VMIN = - VMAX

print("Para - Meta", f_xy(0.31,1.8) )

plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="seismic", levels=500 ) # , norm=Normalize(vmin=VMIN, vmax=VMAX)

np.savetxt(f"{DATA_DIR}/2D_CONTOUR_{EVEC_OUT}_21.dat", f_xy(A0_FINE,WC_FINE))

plt.colorbar(pad=0.01)
#plt.legend()
#plt.xlim(A0_LIST[0],A0_LIST[-1])
#plt.ylim(0)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.savefig(f"{DATA_DIR}/polariton_WCSCAN_A0SCAN_{EVEC_OUT}_21.jpg",dpi=600)
plt.clf()























# Cavity Volume vs. Cavity Frequency

NPTS = 500

f_xy = interp2d(A0_LIST, WC_LIST, (E[2,:,:,0] - E[1,:,:,0]).T, kind='cubic')
WC_FINE = np.linspace(1, 4, NPTS)

#VMIN  = 2*np.pi / (WC_FINE[0]/27.2114) / A0_LIST[-1]**2
#VMIN /= 10**3 / 0.529**3 # nm^3 --> a.u.
#VOLUME_LIST  = np.linspace( VMIN,0.3,NPTS )  # nm**3
VOLUME_LIST  = np.linspace( 0.05,0.3,NPTS )  # nm**3
WV = np.outer( WC_FINE/27.2114, VOLUME_LIST * 10**3 / 0.529**3 )
A0 = np.sqrt( 2 * np.pi / WV )

f_xy_new = np.zeros( (NPTS,NPTS) )
for i in range(NPTS):
    for j in range(NPTS):
        if ( A0[i,j] > 0.4 ): 
            f_xy_new[i,j] = float("Nan")
        else:
            f_xy_new[i,j] = f_xy(A0[i,j],WC_FINE[j])

W,V = np.meshgrid( WC_FINE,  VOLUME_LIST )
plt.contourf( W,V,f_xy_new,levels=1000 )
plt.colorbar(pad=0.01)
plt.xlabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.ylabel("Cavity Volume, $\mathcal{V}$ (nm$^3$)",fontsize=15)
plt.savefig(f"{DATA_DIR}/polariton_VSCAN_WCCSCAN_{EVEC_OUT}_21.jpg",dpi=600)
plt.clf()


np.savetxt(f"{DATA_DIR}/polariton_VSCAN_WCCSCAN_{EVEC_OUT}_21_WGRID.dat", WC_FINE)
np.savetxt(f"{DATA_DIR}/polariton_VSCAN_WCCSCAN_{EVEC_OUT}_21_VGRID.dat", VOLUME_LIST)
np.savetxt(f"{DATA_DIR}/polariton_VSCAN_WCCSCAN_{EVEC_OUT}_21.dat",f_xy_new)






# Cavity Field Strength vs. Cavity Frequency

NPTS = 500

f_xy = interp2d(A0_LIST, WC_LIST, (E[2,:,:,0] - E[1,:,:,0]).T, kind='cubic')
WC_FINE = np.linspace(1, 4, NPTS)
f_xy_new = np.zeros( (NPTS,NPTS) )

FIELD_STRENGTH = np.linspace( 1,20,NPTS )  # V/nm
A0 = np.outer( 1/(WC_FINE/27.2114), FIELD_STRENGTH / 514 ) # A0 = E/WC

for i in range(NPTS):
    for j in range(NPTS):
        if ( A0[i,j] > 0.4 ): 
            f_xy_new[i,j] = float("Nan")
        else:
            f_xy_new[i,j] = f_xy(A0[i,j],WC_FINE[j])

W,F = np.meshgrid( WC_FINE,  FIELD_STRENGTH )
plt.contourf( W,F,f_xy_new.T,levels=1000 )
plt.colorbar(pad=0.01)
plt.xlabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.ylabel("Cavity Field Strength, $|\mathcal{E}|$ (V/nm)",fontsize=15)
plt.savefig(f"{DATA_DIR}/polariton_FSCAN_WCCSCAN_{EVEC_OUT}_21.jpg",dpi=600)
plt.clf()


np.savetxt(f"{DATA_DIR}/polariton_FSCAN_WCCSCAN_{EVEC_OUT}_21_WGRID.dat", WC_FINE)
np.savetxt(f"{DATA_DIR}/polariton_FSCAN_WCCSCAN_{EVEC_OUT}_21_FGRID.dat", FIELD_STRENGTH)
np.savetxt(f"{DATA_DIR}/polariton_FSCAN_WCCSCAN_{EVEC_OUT}_21.dat",f_xy_new)












# Cavity Volume vs. Cavity Frequency

NPTS = 500

f_xy = interp2d(A0_LIST, WC_LIST, (E[0,:,:,0] - E[1,:,:,0]).T, kind='cubic')
WC_FINE = np.linspace(1, 4, NPTS)

#VMIN  = 2*np.pi / (WC_FINE[0]/27.2114) / A0_LIST[-1]**2
#VMIN /= 10**3 / 0.529**3 # nm^3 --> a.u.
#VOLUME_LIST  = np.linspace( VMIN,0.3,NPTS )  # nm**3
VOLUME_LIST  = np.linspace( 0.05,0.3,NPTS )  # nm**3
WV = np.outer( WC_FINE/27.2114, VOLUME_LIST * 10**3 / 0.529**3 )
A0 = np.sqrt( 2 * np.pi / WV )

f_xy_new = np.zeros( (NPTS,NPTS) )
for i in range(NPTS):
    for j in range(NPTS):
        if ( A0[i,j] > 0.4 ): 
            f_xy_new[i,j] = float("Nan")
        else:
            f_xy_new[i,j] = f_xy(A0[i,j],WC_FINE[j])

W,V = np.meshgrid( WC_FINE,  VOLUME_LIST )
plt.contourf( W,V,f_xy_new,levels=1000 )
plt.colorbar(pad=0.01)
plt.xlabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.ylabel("Cavity Volume, $\mathcal{V}$ (nm$^3$)",fontsize=15)
plt.savefig(f"{DATA_DIR}/polariton_VSCAN_WCCSCAN_{EVEC_OUT}_01.jpg",dpi=600)
plt.clf()


np.savetxt(f"{DATA_DIR}/polariton_VSCAN_WCCSCAN_{EVEC_OUT}_01_WGRID.dat", WC_FINE)
np.savetxt(f"{DATA_DIR}/polariton_VSCAN_WCCSCAN_{EVEC_OUT}_01_VGRID.dat", VOLUME_LIST)
np.savetxt(f"{DATA_DIR}/polariton_VSCAN_WCCSCAN_{EVEC_OUT}_01.dat",f_xy_new)






# Cavity Field Strength vs. Cavity Frequency

NPTS = 500

f_xy = interp2d(A0_LIST, WC_LIST, (E[0,:,:,0] - E[1,:,:,0]).T, kind='cubic')
WC_FINE = np.linspace(1, 4, NPTS)
f_xy_new = np.zeros( (NPTS,NPTS) )

FIELD_STRENGTH = np.linspace( 1,20,NPTS )  # V/nm
A0 = np.outer( 1/(WC_FINE/27.2114), FIELD_STRENGTH / 514 ) # A0 = E/WC

for i in range(NPTS):
    for j in range(NPTS):
        if ( A0[i,j] > 0.4 ): 
            f_xy_new[i,j] = float("Nan")
        else:
            f_xy_new[i,j] = f_xy(A0[i,j],WC_FINE[j])

W,F = np.meshgrid( WC_FINE,  FIELD_STRENGTH )
plt.contourf( W,F,f_xy_new.T,levels=1000 )
plt.colorbar(pad=0.01)
plt.xlabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.ylabel("Cavity Field Strength, $|\mathcal{E}|$ (V/nm)",fontsize=15)
plt.savefig(f"{DATA_DIR}/polariton_FSCAN_WCCSCAN_{EVEC_OUT}_01.jpg",dpi=600)
plt.clf()


np.savetxt(f"{DATA_DIR}/polariton_FSCAN_WCCSCAN_{EVEC_OUT}_01_WGRID.dat", WC_FINE)
np.savetxt(f"{DATA_DIR}/polariton_FSCAN_WCCSCAN_{EVEC_OUT}_01_FGRID.dat", FIELD_STRENGTH)
np.savetxt(f"{DATA_DIR}/polariton_FSCAN_WCCSCAN_{EVEC_OUT}_01.dat",f_xy_new)








