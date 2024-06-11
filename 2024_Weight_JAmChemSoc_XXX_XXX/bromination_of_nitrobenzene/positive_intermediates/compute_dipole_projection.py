import numpy as np
from matplotlib import pyplot as plt

NSTATES = 50
DATA_DIR = "MU_PROJECTION"

THETA_1 = 68.8 / 180 * np.pi
PHI_1   = 34.4 / 180 * np.pi
THETA_2 = 68.8 / 180 * np.pi
PHI_2   = 80.2 / 180 * np.pi

E_1_x   = np.sin(THETA_1) * np.cos(PHI_1)
E_1_y   = np.sin(THETA_1) * np.sin(PHI_1)
E_1_z   = np.cos(THETA_1)
E_1     = np.array([E_1_x,E_1_y,E_1_z])

E_2_x   = np.sin(THETA_2) * np.cos(PHI_2)
E_2_y   = np.sin(THETA_2) * np.sin(PHI_2)
E_2_z   = np.cos(THETA_2)
E_2     = np.array([E_2_x,E_2_y,E_2_z])

MU_oBr = np.load("oBr/QCHEM/PLOTS_DATA/DIPOLE_RPA.dat.npy")[:NSTATES,:NSTATES,:]
MU_mBr = np.load("mBr/QCHEM/PLOTS_DATA/DIPOLE_RPA.dat.npy")[:NSTATES,:NSTATES,:]
MU_pBr = np.load("pBr/QCHEM/PLOTS_DATA/DIPOLE_RPA.dat.npy")[:NSTATES,:NSTATES,:]

MU_oBr_2 = np.einsum("abd,bcd->acd", MU_oBr, MU_oBr)
MU_mBr_2 = np.einsum("abd,bcd->acd", MU_mBr, MU_mBr)
MU_pBr_2 = np.einsum("abd,bcd->acd", MU_pBr, MU_pBr)

TMP1   = MU_oBr_2[0,0,:] - MU_mBr_2[0,0,:]
TMP2   = MU_pBr_2[0,0,:] - MU_mBr_2[0,0,:]

print(  TMP1 / np.linalg.norm(TMP1)  )
print(  TMP2 / np.linalg.norm(TMP2)  )

print( E_1 )
print( E_2 )

print( np.dot( TMP1 / np.linalg.norm(TMP1), E_1 ) )
print( np.dot( TMP2 / np.linalg.norm(TMP2), E_2 ) )

exit()

### COMPUTE OVERLAP WITH EACH MOLECULE SEPARATELY ###
# oBr
MU_NORM   = np.zeros( (NSTATES,NSTATES,3) )
MU_PROJ_1 = np.zeros( (NSTATES,NSTATES) )
MU_PROJ_2 = np.zeros( (NSTATES,NSTATES) )
for j in range( NSTATES ):
    for k in range( NSTATES ):
        MU_NORM[j,k,:]   = MU_oBr[j,k,:] / np.linalg.norm( MU_oBr[j,k,:] )
        MU_PROJ_1[j,k]   = np.dot(MU_NORM[j,k,:], E_1 )
        MU_PROJ_2[j,k]   = np.dot(MU_NORM[j,k,:], E_2 )

plt.imshow( np.abs(MU_PROJ_1), origin='lower', cmap="Greys" )
plt.colorbar(pad=0.01)
plt.savefig(f"{DATA_DIR}/MU_oBr_1.jpg",dpi=300)
plt.clf()

plt.imshow( np.abs(MU_PROJ_2), origin='lower', cmap="Greys" )
plt.colorbar(pad=0.01)
plt.savefig(f"{DATA_DIR}/MU_oBr_2.jpg",dpi=300)
plt.clf()


# mBr
MU_NORM   = np.zeros( (NSTATES,NSTATES,3) )
MU_PROJ_1 = np.zeros( (NSTATES,NSTATES) )
MU_PROJ_2 = np.zeros( (NSTATES,NSTATES) )
for j in range( NSTATES ):
    for k in range( NSTATES ):
        MU_NORM[j,k,:]   = MU_mBr[j,k,:] / np.linalg.norm( MU_mBr[j,k,:] )
        MU_PROJ_1[j,k]   = np.dot(MU_NORM[j,k,:], E_1 )
        MU_PROJ_2[j,k]   = np.dot(MU_NORM[j,k,:], E_2 )

plt.imshow( np.abs(MU_PROJ_1), origin='lower', cmap="Greys" )
plt.colorbar(pad=0.01)
plt.savefig(f"{DATA_DIR}/MU_mBr_1.jpg",dpi=300)
plt.clf()

plt.imshow( np.abs(MU_PROJ_2), origin='lower', cmap="Greys" )
plt.colorbar(pad=0.01)
plt.savefig(f"{DATA_DIR}/MU_mBr_2.jpg",dpi=300)
plt.clf()


# pBr
MU_NORM   = np.zeros( (NSTATES,NSTATES,3) )
MU_PROJ_1 = np.zeros( (NSTATES,NSTATES) )
MU_PROJ_2 = np.zeros( (NSTATES,NSTATES) )
for j in range( NSTATES ):
    for k in range( NSTATES ):
        MU_NORM[j,k,:]   = MU_pBr[j,k,:] / np.linalg.norm( MU_pBr[j,k,:] )
        MU_PROJ_1[j,k]   = np.dot(MU_NORM[j,k,:], E_1 )
        MU_PROJ_2[j,k]   = np.dot(MU_NORM[j,k,:], E_2 )

plt.imshow( np.abs(MU_PROJ_1), origin='lower', cmap="Greys" )
plt.colorbar(pad=0.01)
plt.savefig(f"{DATA_DIR}/MU_pBr_1.jpg",dpi=300)
plt.clf()

plt.imshow( np.abs(MU_PROJ_2), origin='lower', cmap="Greys" )
plt.colorbar(pad=0.01)
plt.savefig(f"{DATA_DIR}/MU_pBr_2.jpg",dpi=300)
plt.clf()




### COMPUTE OVERLAP WITH AVERAGE DIPOLE OF EACH MOLECULE SEPARATELY ###
# oBr/mBr
MU_NORM   = np.zeros( (NSTATES,NSTATES,3) )
MU_PROJ_1 = np.zeros( (NSTATES,NSTATES) )
MU_PROJ_2 = np.zeros( (NSTATES,NSTATES) )
for j in range( NSTATES ):
    for k in range( NSTATES ):
        TMP1             = MU_oBr[j,k,:] #/ np.linalg.norm( MU_pBr[j,k,:] )
        TMP2             = MU_mBr[j,k,:] #/ np.linalg.norm( MU_mBr[j,k,:] )
        MU_NORM[j,k,:]   = 0.5 * ( np.abs(TMP1) + np.abs(TMP2) )
        MU_PROJ_1[j,k]   = np.dot(MU_NORM[j,k,:], E_1 )
        MU_PROJ_2[j,k]   = np.dot(MU_NORM[j,k,:], E_2 )

plt.imshow( np.abs(MU_PROJ_1), origin='lower', cmap="Greys" )
plt.colorbar(pad=0.01)
plt.savefig(f"{DATA_DIR}/MU_oBr_mBr_1.jpg",dpi=300)
plt.clf()

plt.imshow( np.abs(MU_PROJ_2), origin='lower', cmap="Greys" )
plt.colorbar(pad=0.01)
plt.savefig(f"{DATA_DIR}/MU_oBr_mBr_2.jpg",dpi=300)
plt.clf()


# pBr/mBr
MU_NORM   = np.zeros( (NSTATES,NSTATES,3) )
MU_PROJ_1 = np.zeros( (NSTATES,NSTATES) )
MU_PROJ_2 = np.zeros( (NSTATES,NSTATES) )
for j in range( NSTATES ):
    for k in range( NSTATES ):
        TMP1             = MU_pBr[j,k,:] #/ np.linalg.norm( MU_pBr[j,k,:] )
        TMP2             = MU_mBr[j,k,:] #/ np.linalg.norm( MU_mBr[j,k,:] )
        MU_NORM[j,k,:]   = 0.5 * ( np.abs(TMP1) + np.abs(TMP2) )
        MU_PROJ_1[j,k]   = np.dot(MU_NORM[j,k,:], E_1 )
        MU_PROJ_2[j,k]   = np.dot(MU_NORM[j,k,:], E_2 )

plt.imshow( np.abs(MU_PROJ_1), origin='lower', cmap="Greys" )
plt.colorbar(pad=0.01)
plt.savefig(f"{DATA_DIR}/MU_pBr_mBr_1.jpg",dpi=300)
plt.clf()

plt.imshow( np.abs(MU_PROJ_2), origin='lower', cmap="Greys" )
plt.colorbar(pad=0.01)
plt.savefig(f"{DATA_DIR}/MU_pBr_mBr_2.jpg",dpi=300)
plt.clf()