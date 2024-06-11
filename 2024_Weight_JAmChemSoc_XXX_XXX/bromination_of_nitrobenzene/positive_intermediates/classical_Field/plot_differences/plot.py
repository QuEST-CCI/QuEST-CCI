import numpy as np
from matplotlib import pyplot as plt
import subprocess as sp

CHI  = np.arange(-15,15,0.5) # V/nm --> a.u.
NCHI = len(CHI)

ENEGRY = np.zeros( (2,NCHI,3) ) # (oBr,mBr), (x,y,z)
E_PF   = np.zeros( (2,NCHI,3) ) # (oBr,mBr), (x,y,z) # wc = 2.8 eV
DIP    = np.zeros( (2,NCHI,3,3) ) # (oBr,mBr), (Ex,Ey,Ez), (dx,dy,dz)

# READ PF CALCULATIONS

for sind,s in enumerate(["oBr","mBr"]):
    for chiind,chi in enumerate( CHI ):
        A0 = chi/514/(2.8/27.2114)
        E_PF[sind,chiind,0] = np.loadtxt(f"../{s}/PF/data_PF/E_THETA_90.0_PHI_0.0_A0_{round(A0,6)}_WC_{round(2.8,6)}_NF_{5}_NM_{50}.dat")[0]
        E_PF[sind,chiind,1] = np.loadtxt(f"../{s}/PF/data_PF/E_THETA_90.0_PHI_90.0_A0_{round(A0,6)}_WC_{round(2.8,6)}_NF_{5}_NM_{50}.dat")[0]
        E_PF[sind,chiind,2] = np.loadtxt(f"../{s}/PF/data_PF/E_THETA_0.0_PHI_0.0_A0_{round(A0,6)}_WC_{round(2.8,6)}_NF_{5}_NM_{50}.dat")[0]
E_PF *= 627.5/27.2114

print( E_PF[0,NCHI//2,0] - E_PF[1,NCHI//2,0] )

# READ FROM GAUSSIAN
for sind,s in enumerate(["oBr","mBr"]):
    for dind,d in enumerate(["X","Y","Z"]):
        for chiind,chi in enumerate( CHI ):
            TMP1 = sp.check_output(f"grep 'SCF Done' ../{s}/DATA_E{d}/STEP_{chiind}/geometry.out | tail -n 1 "+"| awk '{print $5}'", shell=True)
            TMP2 = sp.check_output(f"grep 'Dipole' ../{s}/DATA_E{d}/STEP_{chiind}/geometry.out -A 1 | head -n 2 | tail -n 1 "+"| awk '{print $2, $4, $6}'", shell=True)
            TMP1 = TMP1.decode()
            TMP2 = TMP2.decode().strip().split()
            if ( TMP1 != "" ):
                ENEGRY[sind,chiind,dind] = float( TMP1 )
                DIP[sind,chiind,dind,:]  = np.array(TMP2).astype(float)
            else:
                ENEGRY[sind,chiind,dind] = float("nan")
                DIP[sind,chiind,dind,:]  = float("nan")

# Ortho - Meta
dE_X = ENEGRY[0,:,0] - ENEGRY[1,:,0]
dE_Y = ENEGRY[0,:,1] - ENEGRY[1,:,1]
dE_Z = ENEGRY[0,:,2] - ENEGRY[1,:,2]


# Models
E_0          = ENEGRY[:,NCHI//2,0]
MU_0         = DIP[:,NCHI//2,0,:] / 2.541746
LINEAR_X     = (E_0[0] + CHI[:]/514 * MU_0[0,0]) - (E_0[1] + CHI[:]/514 * MU_0[1,0])
LINEAR_Y     = (E_0[0] + CHI[:]/514 * MU_0[0,1]) - (E_0[1] + CHI[:]/514 * MU_0[1,1])
LINEAR_Z     = (E_0[0] + CHI[:]/514 * MU_0[0,2]) - (E_0[1] + CHI[:]/514 * MU_0[1,2])
QUAD_X       = (E_0[0] + CHI[:]/514 * MU_0[0,0] - 0.5 * (CHI/514)**2 * MU_0[0,0]**2) \
             - (E_0[1] + CHI[:]/514 * MU_0[1,0] - 0.5 * (CHI/514)**2 * MU_0[1,0]**2)
QUAD_Y       = (E_0[0] + CHI[:]/514 * MU_0[0,1] - 0.5 * (CHI/514)**2 * MU_0[0,1]**2) \
             - (E_0[1] + CHI[:]/514 * MU_0[1,1] - 0.5 * (CHI/514)**2 * MU_0[1,1]**2)
QUAD_Z       = (E_0[0] + CHI[:]/514 * MU_0[0,2] - 0.5 * (CHI/514)**2 * MU_0[0,2]**2) \
             - (E_0[1] + CHI[:]/514 * MU_0[1,2] - 0.5 * (CHI/514)**2 * MU_0[1,2]**2)

plt.plot( CHI, E_PF[0,:,0] - E_PF[1,:,0], "o", c='black', alpha=0.25, lw=6, label="E$_x$ (PF)" )
plt.plot( CHI, E_PF[0,:,1] - E_PF[1,:,1], "o", c='red'  , alpha=0.25, lw=6, label="E$_y$ (PF)" )
plt.plot( CHI, E_PF[0,:,2] - E_PF[1,:,2], "o", c='blue' , alpha=0.25, lw=6, label="E$_z$ (PF)" )
plt.plot( CHI, dE_X * 627.5, c='black', lw=3, label="E$_x$ (Gaussian)" )
plt.plot( CHI, dE_Y * 627.5, c='red',   lw=3, label="E$_y$ (Gaussian)" )
plt.plot( CHI, dE_Z * 627.5, c='blue',  lw=3, label="E$_z$ (Gaussian)" )
plt.plot( CHI, LINEAR_X * 627.5, "--", lw=2, c='black', label="E$_x$ (Linear)" )
# plt.plot( CHI, QUAD_X * 627.5, "-.", lw=2, c='black', label="E$_x$ (QUAD)" )
plt.plot( CHI, LINEAR_Y * 627.5, "--", lw=2, c='red', label="E$_y$ (Linear)" )
# plt.plot( CHI, QUAD_Y * 627.5, "-.", lw=2, c='red', label="E$_y$ (QUAD)" )
plt.plot( CHI, LINEAR_Z * 627.5, "--", lw=2, c='blue', label="E$_z$ (Linear)" )
# plt.plot( CHI, QUAD_Z * 627.5, "-.", lw=2, c='blue', label="E$_z$ (QUAD)" )

plt.legend()
plt.xlabel("Electric Field (V/nm)", fontsize=15)
plt.ylabel("$E_\mathrm{ortho} - E_\mathrm{meta}$ (kcal/mol)", fontsize=15)
plt.tight_layout()
plt.savefig("Gaussian_EField.jpg", dpi=300)

print("Done")

