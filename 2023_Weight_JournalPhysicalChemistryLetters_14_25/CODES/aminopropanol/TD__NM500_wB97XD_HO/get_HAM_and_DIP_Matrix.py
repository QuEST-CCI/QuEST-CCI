import numpy as np
import subprocess as sp
from matplotlib import pyplot as plt
import os

NExcStates = 500 # Number of excited states S1, S2, S3, ...


# Need to run Multiwfn ( geometry.fchk 18, 5, geometry.out 2 ) ( 18, 5, geometry.out, 4 )
#sp.call( "python3 ~/Gaussain_scripts/Excited_States/Gen_spectra_gaus16.py" ,shell=True)
sp.call( "./get_dipole_matrix_from_multiwfn.sh", shell=True)



dipole_mat = np.zeros(( NExcStates + 1, NExcStates + 1, 3 )) # All dipole matrix elements
E_adiab = np.zeros( NExcStates + 1 ) # Adiabatic molecular energies

# Read in permanent dipole moments (in a.u.)
permFile = open("dipmom.txt","r").readlines()
dipole_mat[0,0] = np.array( permFile[1].split()[6:9] )
permFile = np.array( permFile[5:NExcStates+5] ) # APPROXIMATE DIPOLES AS HF/CIS DIPOLES
exc_state_perm_dips = np.array( [ permFile[j].split()[1:4] for j in range(NExcStates) ] ).astype(float)
for j in range( NExcStates ):
    dipole_mat[1+j,1+j] = exc_state_perm_dips[j] # CHECK THIS WITH POLAR MOLECULAR SYSTEM

# Read in transition energies (in eV)
Exc_vec = np.array([ permFile[j].split()[4] for j in range(NExcStates) ]).astype(float) # These are in eV

# Get ground state energy from SCF cycle
sp.call(" grep 'SCF Done' geometry.out > GS_Total_Energy.dat ",shell=True) # FOR HF/DFT etc.
#sp.call(" grep 'Wavefunction amplitudes converged' geometry.out > GS_Total_Energy.dat ",shell=True) # FOR CCSD/MP2/etc.
GS_Energy = float(open("GS_Total_Energy.dat","r").readlines()[0].split()[4]) * 27.2114 # Convert to eV

E_adiab[0] = GS_Energy # Shift states by total energy of ground state at this R
for j in range( NExcStates ):
    E_adiab[1+j] = GS_Energy + Exc_vec[j] # In eV

# Read in transition dipole moments
transDipFile = open("transdipmom.txt","r").readlines()
for line in transDipFile:
    t = line.split()
    if (  len(t) == 7 and t[0] != "i" and t[0] != "Transition" ):
        i = int( t[0] )
        j = int( t[1] )
        if ( i <= NExcStates and j <= NExcStates and i != j ): # Exclude terms from ground state and permanant dipoles
            dxyz = np.array( t[2:5] )
            dipole_mat[ i, j ] = dxyz
            dipole_mat[ j, i ] = dxyz

"""
# Read in EOM-CCSD Ground-to-excited dipole moments. These should be exact, but we neglect excited-to-ground transition dipole moments.
sp.call(''' grep "Ground to excited state transition electric dipole moments" geometry.out -A 31 | tail -n 30 > G2E_TdDip.dat ''', shell=True )
transDipFile = np.loadtxt("G2E_TdDip.dat")[:,1:4] # Only <g|r|ej> terms
for j in range(NExcStates):
    dipole_mat[ 0, j+1 ] = transDipFile[j]
    dipole_mat[ j+1, 0 ] = dipole_mat[ 0, j+1 ]
"""

for d in [0,1,2]:
    plt.contourf( np.abs( dipole_mat[:,:,d] ) )
    plt.colorbar()
    plt.title("µ (a.u.)")
    if ( d == 0 ):
        plt.savefig(f"Dipole_Matrix_Ex.jpg")
    if ( d == 1 ):
        plt.savefig(f"Dipole_Matrix_Ey.jpg")
    if ( d == 2 ):
        plt.savefig(f"Dipole_Matrix_Ez.jpg")
    plt.clf()

outDipFile = open(f"dipole_matrix_E.dat","w")
for i in range( NExcStates + 1 ):
    for j in range( i, NExcStates + 1 ):
        outDipFile.write(f"{i} {j} {dipole_mat[i,j,0]} {dipole_mat[i,j,1]} {dipole_mat[i,j,2]}\n")
outDipFile.close()

outExcFile = open("adiabtic_energies.dat","w")
for j in range( NExcStates + 1 ):
    outExcFile.write( f"{j}  {E_adiab[j]}\n" ) # Write in eV




