import numpy as np
from numpy import kron as kron_prod
import sys
import subprocess as sp

##### Braden M. Weight #####
#####  April 17, 2023  #####

# SYNTAX: python3 Pauli-Fierz_DdotE.py A0 WC

def get_globals():
    global NM, NF, A0, wc_eV, wc_AU
    global EVEC_INTS, EVEC_NORM, EVEC_OUT
    global RPA, THETA, PHI

    ##### MAIN USER INPUT SECTION #####
    NM        = 50                  # Number of Electronic States (including ground state)
    NF        = 3                  # Number of Fock Basis States
    RPA       = True                # If True, look for TD-DFT/RPA data rather than TD-DFT/TDA data
    ##### END USER INPUT SECTION  #####
    
    ##### DO NOT MODIFY BELOW HERE #####
    A0    = float(sys.argv[1]) # a.u.
    wc_eV = float(sys.argv[2]) # eV
    
    THETA = 69.4 # theta 2
    PHI   = 79.9 # phi 2

    THETA_RAD = THETA / 180 * np.pi
    PHI_RAD   = PHI / 180 * np.pi

    wc_AU     = wc_eV / 27.2114
    Ex        = np.sin(THETA_RAD) * np.cos(PHI_RAD)
    Ey        = np.sin(THETA_RAD) * np.sin(PHI_RAD)
    Ez        = np.cos(THETA_RAD)
    EVEC = np.array([ Ex,Ey,Ez ]) # Cavity Polarization Vector (input as integers without normalizing)
    EVEC_NORM = EVEC / np.linalg.norm(EVEC)

    sp.call("mkdir -p data_PF", shell=True)

def get_a_op(nf):
    a = np.zeros((nf,nf))
    for m in range(1,nf):
        a[m,m-1] = np.sqrt(m)
    return a.T

def get_H_PF(EAD, MU):
    """
    Input: HAD,   Adiabatic hamiltonian (diagonal) energies from electronic structure
    Output: H_PF, Pauli-Fierz Hamiltonian using Adiabatic-Fock basis states
    """

    print (f"Dimension = {(NM*NF)}")

    I_ph = np.identity(NF)
    I_el = np.identity(NM)
    a_op = get_a_op(NF)
    q_op = a_op.T + a_op
    MU   = np.einsum("d,JKd->JK", EVEC_NORM[:], MU[:,:,:] )

    H_EL = np.diag( EAD )
    H_PH = np.diag( np.arange(NF) * wc_AU )

    H_PF   = kron_prod(H_EL, I_ph)                        # Matter
    H_PF  += kron_prod(I_el, H_PH)                        # Photon

    H_PF  += kron_prod(MU, q_op) * wc_AU * A0             # Interaction
    H_PF  += kron_prod(MU @ MU, I_ph) * wc_AU * A0**2     # Dipole Self-Energy

    if ( H_PF.size * H_PF.itemsize * 10 ** -9 > 0.5 ): # If H_PF larger than 0.5 GB, tell the user.
        print(f"\tWARNING!!! Large Matrix" )
        print(f"\tMemory size of numpy array in (MB, GB): ({round(H_PF.size * H_PF.itemsize * 10 ** -6,2)},{round(H_PF.size * H_PF.itemsize * 10 ** -9,2)})" )

    return H_PF

def get_ADIABATIC_DATA():

    EAD  = np.zeros(( NM ))
    MU  = np.zeros(( NM,NM,3 ))

    if ( RPA == True ):
        EAD += np.loadtxt(f"../../TDDFT_1000/PLOTS_DATA/ADIABATIC_ENERGIES_RPA.dat")[:NM] # in AU
        MU  += np.load(f"../../TDDFT_1000/PLOTS_DATA/DIPOLE_RPA.dat.npy")[:NM,:NM] # in AU
    else:
        EAD += np.loadtxt(f"../../TDDFT_1000/PLOTS_DATA/ADIABATIC_ENERGIES_TDA.dat")[:NM] # in AU
        MU  += np.load(f"../../TDDFT_1000/PLOTS_DATA/DIPOLE_TDA.dat.npy")[:NM,:NM] # in AU

    return EAD, MU

def SolvePlotandSave(H_PF,EAD,MU):

        # Diagonalize Hamiltonian
        E, U = np.linalg.eigh( H_PF ) # This is exact solution --> Should we ever try Davidson diagonalization ?
        
        # Save Data
        np.savetxt( f"data_PF/E_THETA_{round(THETA,2)}_PHI_{round(PHI,2)}_A0_{round(A0,6)}_WC_{round(wc_eV,6)}_NF_{NF}_NM_{NM}.dat", E * 27.2114 )
        #np.savetxt( f"data_PF/U_THETA_{round(THETA,2)}_PHI_{round(PHI,2)}_A0_{round(A0,6)}_WC_{round(wc_eV,6)}_NF_{NF}_NM_{NM}.dat", U ) # These can be large
        np.save( f"data_PF/U_THETA_{round(THETA,2)}_PHI_{round(PHI,2)}_A0_{round(A0,6)}_WC_{round(wc_eV,6)}_NF_{NF}_NM_{NM}.dat", U ) # Binary is smaller

        print ( A0, wc_eV )

        # Save original EAD and MU for convenience
        np.savetxt( f"data_PF/EAD.dat", EAD * 27.2114 )
        np.save( f"data_PF/MU.dat", MU )

def main():

    get_globals()

    EAD, MU = get_ADIABATIC_DATA() 
    H_PF    = get_H_PF( EAD, MU )
    SolvePlotandSave( H_PF, EAD, MU)

    


if __name__ == "__main__":
    main()