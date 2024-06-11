import numpy as np
from matplotlib import pyplot as plt

A0_LIST = np.array([0.01,0.1,0.2,0.3,0.4])

NF_FIX  = 5
NM_LIST = np.array([2,3,4,5,10,25,50,75,100,200,300,400,500,600,700,800,900,1000])
NM_FIX  = 50
NF_LIST = np.array([2,3,4,5,6,7,8,9,10,20,50])

E_MSCAN = np.zeros( (len(A0_LIST),len(NM_LIST)) )
E_FSCAN = np.zeros( (len(A0_LIST),len(NF_LIST)) )

for A0i,A0 in enumerate(A0_LIST):
    for NMi,NM in enumerate(NM_LIST):
        E_MSCAN[A0i,NMi] = np.loadtxt( f"PF_NF{NF_FIX}_NM{NM}/data_PF/E_0_1_0_A0_{A0}_WC_5.0_NF_{NF_FIX}_NM_{NM}.dat" )[0] / 27.2114 * 627.5

    for NFi,NF in enumerate(NF_LIST):
        E_FSCAN[A0i,NFi] = np.loadtxt( f"PF_NF{NF}_NM{NM_FIX}/data_PF/E_0_1_0_A0_{A0}_WC_5.0_NF_{NF}_NM_{NM_FIX}.dat" )[0] / 27.2114 * 627.5

### MATTER SCAN ###
for A0i,A0 in enumerate(A0_LIST):
    plt.plot( NM_LIST, E_MSCAN[A0i] - np.min(E_MSCAN[A0i,0]), "-o", label="$A_0$"+f" = {A0}" )
plt.legend()
plt.xlabel("Number of Electronic States",fontsize=15)
plt.ylabel("Energy (kcal/mol)",fontsize=15)
plt.title("$N_\mathrm{F} = $"+f"{NF_FIX}",fontsize=15)
plt.savefig(f"NM_SCAN_NF_{NF_FIX}.jpg",dpi=300)
plt.clf()

### FOCK SCAN ###
for A0i,A0 in enumerate(A0_LIST):
    plt.semilogx( NF_LIST, E_FSCAN[A0i] - np.min(E_FSCAN[A0i,0]), "-o", label="$A_0$"+f" = {A0}" )
plt.legend()
plt.xlabel("Number of Fock States",fontsize=15)
plt.ylabel("Energy (kcal/mol)",fontsize=15)
plt.title("$N_\mathrm{el} = $"+f"{NM_FIX}",fontsize=15)
plt.savefig(f"NF_SCAN_NM_{NM_FIX}.jpg",dpi=300)
plt.clf()



