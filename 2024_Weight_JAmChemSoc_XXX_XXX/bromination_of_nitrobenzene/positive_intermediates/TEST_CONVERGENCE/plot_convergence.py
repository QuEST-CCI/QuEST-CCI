import numpy as np
from matplotlib import pyplot as plt

A0_LIST = np.array([0.01,0.1,0.2,0.3])

NF_FIX  = 10
#NM_LIST = np.array([2,3,4,5,10,25,50,75,100,200,300,400,500,600,700,800,900,1000])
NM_LIST = np.array([2,3,4,5,10,25,50,75,100,200,300])
NM_FIX  = 100
NF_LIST = np.array([2,3,4,5,6,7,8,9,10,20,50])
    
THETA = 69.4 # theta 2
PHI   = 79.9 # phi 2

E_MSCAN = np.zeros( (2, len(A0_LIST),len(NM_LIST)) )
E_FSCAN = np.zeros( (2, len(A0_LIST),len(NF_LIST)) )


# READ IN mBr
for A0i,A0 in enumerate(A0_LIST):
    print("Reading in mBr...")
    for NMi,NM in enumerate(NM_LIST):
        E_MSCAN[0,A0i,NMi] = np.loadtxt( f"../mBr/QCHEM/TEST_CONVERGENCE/NEW_CONVERGENCE/PF_NF{NF_FIX}_NM{NM}/data_PF/E_THETA_{round(THETA,2)}_PHI_{round(PHI,2)}_A0_{A0}_WC_1.8_NF_{NF_FIX}_NM_{NM}.dat" )[0] / 27.2114 * 627.5

    for NFi,NF in enumerate(NF_LIST):
        E_FSCAN[0,A0i,NFi] = np.loadtxt( f"../mBr/QCHEM/TEST_CONVERGENCE/NEW_CONVERGENCE/PF_NF{NF}_NM{NM_FIX}/data_PF/E_THETA_{round(THETA,2)}_PHI_{round(PHI,2)}_A0_{A0}_WC_1.8_NF_{NF}_NM_{NM_FIX}.dat" )[0] / 27.2114 * 627.5

# READ IN pBr
for A0i,A0 in enumerate(A0_LIST):
    print("Reading in pBr...")
    for NMi,NM in enumerate(NM_LIST):
        E_MSCAN[1,A0i,NMi] = np.loadtxt( f"../pBr/QCHEM/TEST_CONVERGENCE/NEW_CONVERGENCE/PF_NF{NF_FIX}_NM{NM}/data_PF/E_THETA_{round(THETA,2)}_PHI_{round(PHI,2)}_A0_{A0}_WC_1.8_NF_{NF_FIX}_NM_{NM}.dat" )[0] / 27.2114 * 627.5

    for NFi,NF in enumerate(NF_LIST):
        E_FSCAN[1,A0i,NFi] = np.loadtxt( f"../pBr/QCHEM/TEST_CONVERGENCE/NEW_CONVERGENCE/PF_NF{NF}_NM{NM_FIX}/data_PF/E_THETA_{round(THETA,2)}_PHI_{round(PHI,2)}_A0_{A0}_WC_1.8_NF_{NF}_NM_{NM_FIX}.dat" )[0] / 27.2114 * 627.5


dETS_MSCAN = E_MSCAN[1,:,:] - E_MSCAN[0,:,:] # Para - Meta
dETS_FSCAN = E_FSCAN[1,:,:] - E_FSCAN[0,:,:] # Para - Meta


### MATTER SCAN ###
for A0i,A0 in enumerate(A0_LIST):
    #plt.semilogx( NM_LIST, dETS_MSCAN[A0i] - np.min(dETS_MSCAN[A0i,0]), "-o", label="$A_0$"+f" = {A0}" )
    plt.semilogx( NM_LIST, dETS_MSCAN[A0i] - np.min(dETS_MSCAN[A0i,-1]), "-o", label="$A_0$"+f" = {A0}" )
    #plt.plot( NM_LIST, dETS_MSCAN[A0i], "-o", label="$A_0$"+f" = {A0}" )
plt.legend()
plt.xlabel("Number of Electronic States",fontsize=15)
plt.ylabel("Energy (kcal/mol)",fontsize=15)
plt.title("$N_\mathrm{F} = $"+f"{NF_FIX}",fontsize=15)
plt.savefig(f"NM_SCAN_NF_{NF_FIX}.jpg",dpi=300)
plt.clf()

### FOCK SCAN ###
for A0i,A0 in enumerate(A0_LIST):
    plt.plot( NF_LIST, dETS_FSCAN[A0i] - np.min(dETS_FSCAN[A0i,-1]), "-o", label="$A_0$"+f" = {A0}" )
    #plt.plot( NF_LIST, dETS_FSCAN[A0i], "-o", label="$A_0$"+f" = {A0}" )
plt.legend()
plt.xlabel("Number of Fock States",fontsize=15)
plt.ylabel("Energy (kcal/mol)",fontsize=15)
plt.title("$N_\mathrm{el} = $"+f"{NM_FIX}",fontsize=15)
plt.savefig(f"NF_SCAN_NM_{NM_FIX}.jpg",dpi=300)
plt.clf()



