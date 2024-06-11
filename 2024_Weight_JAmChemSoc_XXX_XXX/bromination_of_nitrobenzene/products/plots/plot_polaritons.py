import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp2d
import os

A0_LIST = np.arange( 0, 0.4+0.01, 0.01 )
WC_LIST = np.arange( 0, 20.0+1.0, 1.0 )

NM = 50
NF = 5

EVEC_INTS = np.array([1,0,0])
EVEC_NORM = EVEC_INTS / np.linalg.norm(EVEC_INTS)
EVEC_OUT = "_".join(map(str,EVEC_INTS))

NA0  = len(A0_LIST)
NWC  = len(WC_LIST)
NPOL = NM * NF

DIRS = ["oBr","mBr","pBr"]
EPOL = np.zeros(( len(DIRS), NPOL, NA0, NWC )) # [oBr,mBr]
PHOT_0 = np.zeros(( len(DIRS), NA0, NWC )) # (o,m,p) x NPOL
EXCI_0 = np.zeros(( len(DIRS), NA0, NWC )) # (o,m,p) x NPOL
#U = np.zeros(( len(DIRS), NA0, NWC, NPOL, NPOL )) # (o,m,p) x NPOL
MU       = np.zeros(( len(DIRS), NM, NM ))
MU_POL_0 = np.zeros(( len(DIRS), NA0, NWC ))

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

I_m  = np.identity(NM)
I_ph = np.identity(NF)
PHOT_OP_FULL = np.kron( I_m, N_a_OP )
EXCI_OP_FULL = np.kron( N_c_OP, I_ph  )

if ( not os.path.isfile("./E_POL.npy") or \
     not os.path.isfile("./U_POL.npy") or \
     not os.path.isfile("./PHOT_0_POL.npy") or \
     not os.path.isfile("./MU.npy") or \
     not os.path.isfile("./MU_POL_0.npy") or \
     not os.path.isfile("./EXCI_0_POL.npy") ):

    for GEOMind,GEOM in enumerate(DIRS):
        TMP = np.load(f"../{GEOM}/QCHEM/PLOTS_DATA/DIPOLE_RPA.dat.npy")[:NM,:NM,:]
        MU[GEOMind,:,:] = np.einsum("d,JKd->JK", EVEC_INTS, TMP)

    MU_POL_OP = np.array([  np.kron( MU[0,:,:], I_ph ),\
                            np.kron( MU[1,:,:], I_ph ),\
                            np.kron( MU[2,:,:], I_ph ) ])

    for A0_IND,A0 in enumerate( A0_LIST ):
        for WC_IND,WC in enumerate( WC_LIST ):
            A0 = round(A0,4)
            WC = round(WC,4)

            for GEOMind,GEOM in enumerate(DIRS):

                EPOL[GEOMind,:,A0_IND,WC_IND] = np.loadtxt(f"../{GEOM}/QCHEM/PF_NF{NF}_NM{NM}/data_PF/E_{EVEC_OUT}_A0_{A0}_WC_{WC}_NF_{NF}_NM_{NM}.dat") / 27.2114 * 630
                U0 = np.load(f"../{GEOM}/QCHEM/PF_NF{NF}_NM{NM}/data_PF/U_{EVEC_OUT}_A0_{A0}_WC_{WC}_NF_{NF}_NM_{NM}.dat.npy")[:,0]
                PHOT_0[GEOMind,A0_IND,WC_IND]   = U0 @ PHOT_OP_FULL[:,:] @ U0
                EXCI_0[GEOMind,A0_IND,WC_IND]   = U0 @ EXCI_OP_FULL[:,:] @ U0
                MU_POL_0[GEOMind,A0_IND,WC_IND] = U0 @ MU_POL_OP[GEOMind,:,:]  @ U0


            print(A0, WC, np.round(PHOT_0[0,A0_IND,WC_IND],5), round(EXCI_0[0,A0_IND,WC_IND],5), np.round(PHOT_0[1,A0_IND,WC_IND],5), round(EXCI_0[1,A0_IND,WC_IND],5) )
            print( "\tDIP:", np.round(MU_POL_0[0,A0_IND,WC_IND],5), np.round(MU[0,0,0],5), np.round(MU_POL_0[1,A0_IND,WC_IND],5), np.round(MU[1,0,0],5) )
    

    np.save("E_POL.npy", EPOL)
    np.save("MU_POL_0.npy", MU_POL_0)
    #np.save("U_POL.npy", U)
    np.save("PHOT_0_POL.npy", PHOT_0)
    np.save("EXCI_0_POL.npy", EXCI_0)

else:
    print("Found data. Reading it.")
    EPOL      = np.load("E_POL.npy")
    MU        = np.load("MU.npy")
    MU_POL_0  = np.load("MU_POL_0.npy")
    #U      = np.load("U_POL.npy")
    PHOT_0    = np.load("PHOT_0_POL.npy")
    EXCI_0    = np.load("EXCI_0_POL.npy")


# SAVE GROUND STATE ENERGIES FOR ALL PARAMETERS TO FILE
print(EPOL.shape) # EPOL.shape = (2, NPOL, NA0, NWC)
np.savetxt("E0_ORTHO.dat", EPOL[0,0,:,:] ) # meV --> Kcal/mol
np.savetxt("E0_META.dat",  EPOL[1,0,:,:] ) # meV --> Kcal/mol
np.savetxt("E0_ORTHO_EZERO.dat", (EPOL[0,0,:,:] - np.min( EPOL[1,:] )) ) # meV --> Kcal/mol
np.savetxt("E0_META_EZERO.dat",  (EPOL[1,0,:,:] - np.min( EPOL[1,:] )) ) # meV --> Kcal/mol
np.savetxt("dE0_ORTHO-META.dat", (EPOL[0,0,:,:] - EPOL[1,0,:,:]) ) # meV --> Kcal/mol
#exit()



######## PLOT 5 (2D) DIP oBr ########

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, np.abs(MU_POL_0[0,:,:].T), kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1], 500)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1], 500)

X,Y = np.meshgrid(A0_FINE,WC_FINE)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="hot", levels=500)

np.savetxt("2D_CONTOUR_DIP_oBr.dat", f_xy(A0_LIST,WC_LIST))

plt.colorbar(pad=0.01)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.title("oBr Ground State Dipole, $|\mu_{00}| (a.u.)$",fontsize=15)
plt.savefig(f"DIP_wcSCAN_A0SCAN_oBr.jpg",dpi=600)
plt.clf()


######## PLOT 5 (2D) DIP mBr ########

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, np.abs(MU_POL_0[1,:,:].T), kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1], 500)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1], 500)

X,Y = np.meshgrid(A0_FINE,WC_FINE)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="hot", levels=500)

np.savetxt("2D_CONTOUR_DIP_mBr.dat", f_xy(A0_LIST,WC_LIST))

plt.colorbar(pad=0.01)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.title("mBr Ground State Dipole, $|\mu_{00}|$ (a.u.)",fontsize=15)
plt.savefig(f"DIP_wcSCAN_A0SCAN_mBr.jpg",dpi=600)
plt.clf()




######## PLOT 5 (2D) PHOT oBr ########

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, PHOT_0[0,:,:].T, kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1], 500)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1], 500)

X,Y = np.meshgrid(A0_FINE,WC_FINE)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="hot_r", levels=1000)

np.savetxt("2D_CONTOUR_PHOT_oBr.dat", f_xy(A0_LIST,WC_LIST))

plt.colorbar(pad=0.01)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.title("oBr Average Photon Number, $ \langle \hat{a}^\dag \hat{a} \\rangle $",fontsize=15)
plt.savefig(f"PHOT_wcSCAN_A0SCAN_oBr.jpg",dpi=600)
plt.clf()


######## PLOT 5 (2D) PHOT mBr ########

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, PHOT_0[1,:,:].T, kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1], 500)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1], 500)

X,Y = np.meshgrid(A0_FINE,WC_FINE)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="hot_r", levels=1000)

np.savetxt("2D_CONTOUR_PHOT_mBr.dat", f_xy(A0_LIST,WC_LIST))

plt.colorbar(pad=0.01)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.title("mBr Average Photon Number, $ \langle  \hat{a}^\dag \hat{a} \\rangle $",fontsize=15)
plt.savefig(f"PHOT_wcSCAN_A0SCAN_mBr.jpg",dpi=600)
plt.clf()


######## PLOT 5 (2D) EXCI oBr ########

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, EXCI_0[0,:,:].T, kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1], 500)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1], 500)

X,Y = np.meshgrid(A0_FINE,WC_FINE)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="hot_r", levels=1000)

np.savetxt("2D_CONTOUR_EXCI_oBr.dat", f_xy(A0_LIST,WC_LIST))

plt.colorbar(pad=0.01)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.title("oBr Average Exciton Number, $ \langle \hat{c}^\dag \hat{c} \\rangle $",fontsize=15)
plt.savefig(f"EXCI_wcSCAN_A0SCAN_oBr.jpg",dpi=600)
plt.clf()


######## PLOT 5 (2D) EXCI mBr ########

# Let's interpolate the coarse 2D grid data
f_xy = interp2d(A0_LIST, WC_LIST, EXCI_0[1,:,:].T, kind='cubic')
# Define fine grid
A0_FINE = np.linspace(A0_LIST[0], A0_LIST[-1], 500)
WC_FINE = np.linspace(WC_LIST[0], WC_LIST[-1], 500)

X,Y = np.meshgrid(A0_FINE,WC_FINE)
plt.contourf( X,Y, f_xy(A0_FINE,WC_FINE), cmap="hot_r", levels=1000)

np.savetxt("2D_CONTOUR_EXCI_mBr.dat", f_xy(A0_LIST,WC_LIST))

plt.colorbar(pad=0.01)
plt.xlabel("Coupling Strength, A$_0$ (a.u.)",fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_c$ (eV)",fontsize=15)
plt.title("mBr Average Exciton Number, $ \langle \hat{c}^\dag \hat{c} \\rangle$",fontsize=15)
plt.savefig(f"EXCI_wcSCAN_A0SCAN_mBr.jpg",dpi=600)
plt.clf()




















EZERO  = EPOL[2,0,0,0]

# E as functions of A0
colors = plt.cm.brg(np.linspace(0, 1, NWC))
for WC_IND,WC in enumerate( WC_LIST ):
    WC    = round(WC,4)
    if ( WC_IND in [0,NWC-1] ):
        plt.plot( A0_LIST, EPOL[0,0,:,WC_IND] - EZERO, "-", lw=4, c=colors[WC_IND], label=f"WC = {WC} (oBr)" )
        plt.plot( A0_LIST, EPOL[1,0,:,WC_IND] - EZERO, ".", lw=4, c=colors[WC_IND], label=f"WC = {WC} (mBr)" )
        plt.plot( A0_LIST, EPOL[2,0,:,WC_IND] - EZERO, "o", lw=4, c=colors[WC_IND], label=f"WC = {WC} (pBr)" )
    else:
        plt.plot( A0_LIST, EPOL[0,0,:,WC_IND] - EZERO, "-", lw=4, c=colors[WC_IND] )
        plt.plot( A0_LIST, EPOL[1,0,:,WC_IND] - EZERO, ".", c=colors[WC_IND] )
        plt.plot( A0_LIST, EPOL[2,0,:,WC_IND] - EZERO, "o", c=colors[WC_IND] )

plt.xlabel("Light-Matter Coupling Strength, $A_0$ (a.u.)", fontsize=15)
plt.ylabel("Polaritonic Ground State Energy (kcal/mol)", fontsize=15)
plt.xlim(A0_LIST[0],A0_LIST[-1])
plt.legend()
plt.savefig("E0_A0SCAN.jpg", dpi=600)
plt.clf()


# E as functions of WC
colors = plt.cm.brg(np.linspace(0, 1, NA0))
for A0_IND,A0 in enumerate( A0_LIST ):
    A0 = round(A0,4)
    if ( A0_IND in [0,NA0-1] ):
        plt.plot( WC_LIST, EPOL[0,0,A0_IND,:] - EZERO, "-", lw=4, c=colors[A0_IND], label=f"A0 = {A0} (oBr)" )
        plt.plot( WC_LIST, EPOL[1,0,A0_IND,:] - EZERO, ".", lw=4, c=colors[A0_IND], label=f"A0 = {A0} (mBr)" )
        plt.plot( WC_LIST, EPOL[2,0,A0_IND,:] - EZERO, "o", lw=4, c=colors[A0_IND], label=f"A0 = {A0} (pBr)" )
    else:
        plt.plot( WC_LIST, EPOL[0,0,A0_IND,:] - EZERO, "-", lw=4, c=colors[A0_IND] )
        plt.plot( WC_LIST, EPOL[1,0,A0_IND,:] - EZERO, ".", c=colors[A0_IND] )
        plt.plot( WC_LIST, EPOL[2,0,A0_IND,:] - EZERO, "o", c=colors[A0_IND] )

plt.xlabel("Cavity Frequency, $\omega_C$ (eV)", fontsize=15)
plt.ylabel("Polaritonic Ground State Energy (meV)", fontsize=15)
plt.xlim(WC_LIST[0],WC_LIST[-1])
plt.legend()
plt.savefig("E0_WCSCAN.jpg", dpi=600)
plt.clf()




# dE as functions of A0
colors = plt.cm.brg(np.linspace(0, 1, NWC))
for WC_IND,WC in enumerate( WC_LIST ):
    WC    = round(WC,4)
    if ( WC_IND in [0,NWC-1] ):
        plt.plot( A0_LIST, EPOL[0,0,:,WC_IND] - EPOL[2,0,:,WC_IND], "-", lw=4, c=colors[WC_IND], label=f"WC = {WC} eV (o-p)" )
        plt.plot( A0_LIST, EPOL[0,0,:,WC_IND] - EPOL[2,0,:,WC_IND], "o", lw=4, c=colors[WC_IND], label=f"WC = {WC} eV (m-p)" )
    else:
        plt.plot( A0_LIST, EPOL[0,0,:,WC_IND] - EPOL[2,0,:,WC_IND], "-", lw=4, c=colors[WC_IND] )
        plt.plot( A0_LIST, EPOL[1,0,:,WC_IND] - EPOL[2,0,:,WC_IND], "o", lw=4, c=colors[WC_IND] )


plt.xlabel("Light-Matter Coupling Strength, $A_0$ (a.u.)", fontsize=15)
plt.ylabel("dE(o-p); dE(m-p) (kcal/mol)", fontsize=15)
plt.xlim(A0_LIST[0],A0_LIST[-1])
plt.legend()
plt.savefig("dE0_A0SCAN_02_12.jpg", dpi=600)
plt.clf()


# dE as functions of WC
colors = plt.cm.brg(np.linspace(0, 1, NA0))
for A0_IND,A0 in enumerate( A0_LIST ):
    A0 = round(A0,4)
    if ( A0_IND in [0,NA0-1] ):
        plt.plot( WC_LIST, EPOL[0,0,A0_IND,:] - EPOL[2,0,A0_IND,:], "-", lw=4, c=colors[A0_IND], label=f"A0 = {A0} a.u. (o-p)" )
        plt.plot( WC_LIST, EPOL[1,0,A0_IND,:] - EPOL[2,0,A0_IND,:], "o", lw=4, c=colors[A0_IND], label=f"A0 = {A0} a.u. (m-p)" )
    else:
        plt.plot( WC_LIST, EPOL[0,0,A0_IND,:] - EPOL[2,0,A0_IND,:], "-", lw=4, c=colors[A0_IND] )
        plt.plot( WC_LIST, EPOL[1,0,A0_IND,:] - EPOL[2,0,A0_IND,:], "o", lw=4, c=colors[A0_IND] )

plt.xlabel("Cavity Frequency, $\omega_C$ (eV)", fontsize=15)
plt.ylabel("dE(o-p); dE(m-p) (kcal/mol)", fontsize=15)
plt.xlim(WC_LIST[0],WC_LIST[-1])
plt.legend()
plt.savefig("dE0_WCSCAN_02_12.jpg", dpi=600)
plt.clf()


# A0/WC SCAN 2D CONTOUR
plt.contourf( A0_LIST, WC_LIST, EPOL[0,0,:,:].T - EPOL[2,0,:,:].T, cmap="seismic", levels=1000 )
plt.xlabel("Light-Matter Coupling Strength, $A_0$ (a.u.)", fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_C$ (eV)", fontsize=15)
plt.title("dE(o-p), (kcal/mol)", fontsize=15)
plt.colorbar(pad=0.01)
plt.savefig("dE_A0_WC_02.jpg", dpi=600)
plt.clf()

# A0/WC SCAN 2D CONTOUR
plt.contourf( A0_LIST, WC_LIST, EPOL[1,0,:,:].T - EPOL[2,0,:,:].T, cmap="seismic", levels=1000 )
plt.xlabel("Light-Matter Coupling Strength, $A_0$ (a.u.)", fontsize=15)
plt.ylabel("Cavity Frequency, $\omega_C$ (eV)", fontsize=15)
plt.title("dE(m-p), (kcal/mol)", fontsize=15)
plt.colorbar(pad=0.01)
plt.savefig("dE_A0_WC_12.jpg", dpi=600)
plt.clf()
