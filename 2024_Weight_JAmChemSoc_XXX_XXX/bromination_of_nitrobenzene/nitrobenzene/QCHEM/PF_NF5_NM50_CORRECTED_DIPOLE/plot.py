import numpy as np
from matplotlib import pyplot as plt

A0_LIST = np.arange( 0.0, 0.31, 0.01 )
WC_LIST = np.array([1.8])

NPOL = 50*5

ENERGY = np.zeros( (2, len(A0_LIST), len(WC_LIST)) ) # 100, 010
WFN    = np.zeros( (len(A0_LIST), len(WC_LIST), NPOL) ) # 100, 010
MU_X   = np.load("../PLOTS_DATA/DIPOLE_RPA.dat.npy")[:50,:50,0]
MU_Y   = np.load("../PLOTS_DATA/DIPOLE_RPA.dat.npy")[:50,:50,1]

for A0i,A0 in enumerate( A0_LIST ):
    for WCi,WC in enumerate( WC_LIST ):
        A0 = round(A0,4)
        WC = round(WC,4)
        ENERGY[0,A0i,WCi] = np.loadtxt( f"data_PF/E_1_0_0_A0_{A0}_WC_{WC}_NF_5_NM_50.dat" )[0] / 27.2114 * 630
        ENERGY[1,A0i,WCi] = np.loadtxt( f"data_PF/E_0_1_0_A0_{A0}_WC_{WC}_NF_5_NM_50.dat" )[0] / 27.2114 * 630

        WFN[A0i,WCi,:] = np.load( f"data_PF/U_1_0_0_A0_{A0}_WC_{WC}_NF_5_NM_50.dat.npy" )[:,0]
        WFN[A0i,WCi,:] = np.load( f"data_PF/U_0_1_0_A0_{A0}_WC_{WC}_NF_5_NM_50.dat.npy" )[:,0]

for WCi,WC in enumerate( WC_LIST ):
    WC = round(WC,4)
    plt.plot( A0_LIST, ENERGY[0,:,WCi] - ENERGY[0,0,WCi], "-", c="black", lw=4, label="X" )
    plt.plot( A0_LIST, ENERGY[1,:,WCi] - ENERGY[1,0,WCi], "-", c="red", lw=4, label="Y" )

    plt.legend()
    plt.xlim(A0_LIST[0], A0_LIST[-1])
    plt.xlabel("Coupling Strength, $A_0$ (a.u.)", fontsize=15)
    plt.ylabel("Energy (kcal/mol)", fontsize=15)
    plt.tight_layout()
    plt.savefig(f"E_GS_WC_{WC}.jpg", dpi=300)
    plt.clf()


for WCi,WC in enumerate( WC_LIST ):
    WC = round(WC,4)
    plt.plot( A0_LIST, ENERGY[0,:,WCi] - ENERGY[1,:,WCi], "-", c="black", lw=4, label="X-Y" )

    plt.legend()
    plt.xlabel("Coupling Strength, $A_0$ (a.u.)", fontsize=15)
    plt.ylabel("Energy (kcal/mol)", fontsize=15)
    plt.xlim(A0_LIST[0], A0_LIST[-1])
    plt.tight_layout()
    plt.savefig(f"dE_GS_WC_{WC}.jpg", dpi=300)
    plt.clf()



# GS Dipole

MU_X_BARE_FULL = np.kron( MU_X[:,:], np.identity(5) )
for WCi,WC in enumerate( WC_LIST ):
    WC = round(WC,4)

    MU_X_POL = np.einsum( "as,sk,ak->a", WFN[:,WCi,:], MU_X_BARE_FULL[:,:], WFN[:,WCi,:] )
    plt.plot( A0_LIST, A0_LIST*0 + MU_X[0,0], "--", c="black", lw=2 )
    plt.plot( A0_LIST, MU_X_POL, "-", c="black", lw=4, label="X" )

    plt.legend()
    plt.xlabel("Coupling Strength, $A_0$ (a.u.)", fontsize=15)
    plt.ylabel("Dipole Moment, $\mu_{00}$ (a.u.)", fontsize=15)
    plt.xlim(A0_LIST[0], A0_LIST[-1])
    plt.tight_layout()
    plt.savefig(f"DIPOLE_GS_WC_{WC}_X.jpg", dpi=300)
    plt.clf()


MU_Y_BARE_FULL = np.kron( MU_Y[:,:], np.identity(5) )
for WCi,WC in enumerate( WC_LIST ):
    WC = round(WC,4)

    MU_Y_POL = np.einsum( "as,sk,ak->a", WFN[:,WCi,:], MU_Y_BARE_FULL[:,:], WFN[:,WCi,:] )
    plt.plot( A0_LIST, A0_LIST*0 + MU_Y[0,0], "--", c="red", lw=2 )
    plt.plot( A0_LIST, MU_Y_POL, "-", c="red", lw=4, label="Y" )

    plt.legend()
    plt.xlabel("Coupling Strength, $A_0$ (a.u.)", fontsize=15)
    plt.ylabel("Dipole Moment, $\mu_{00}$ (a.u.)", fontsize=15)
    plt.xlim(A0_LIST[0], A0_LIST[-1])
    plt.tight_layout()
    plt.savefig(f"DIPOLE_GS_WC_{WC}_Y.jpg", dpi=300)
    plt.clf()




