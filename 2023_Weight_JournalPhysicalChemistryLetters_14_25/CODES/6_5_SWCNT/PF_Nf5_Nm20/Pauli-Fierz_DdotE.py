import numpy as np
from scipy.linalg import eigh as SC_EXACT_EIGH
from numpy import kron as ꕕ
import sys
from matplotlib import pyplot as plt
import multiprocessing as mp

def ĉ(nf):
    a = np.zeros((nf,nf))
    for m in range(1,nf):
        a[m,m-1] = np.sqrt(m)
    return a.T

class params:

    Nad = 20 # Number of Electronic Basis Excited States
    nf = 5 # Number of Fock Basis States

    η_1 = float( sys.argv[1] )
    ωc  = 1.0    # eV

    χ_1 = (ωc/27.2114) * η_1

def getĤpol(Had, µ, params):
    """
    Input: Had, adiabatic hamiltonian (diagonal) energies from electronic structure
    Output: Hpol, Hamiltonian after Kron. product with single- or double-mode photon cavity space
    """

    ns = len(Had)
    nf = params.nf
    ωc = params.ωc 
    χ_1 = params.χ_1

    print (f"Dimension = {(ns*nf)}")

    #------------------------
    Iₚ_1 = np.identity(nf)
    Iₘ = np.identity(ns)
    #------ Photonic Part -----------------------
    Hₚ_1 = np.identity(nf)
    Hₚ_1[np.diag_indices(nf)] = np.arange(nf) * (ωc/27.2114)
    â = ĉ(nf) 
    #--------------------------------------------
    #       matter ⊗ photon 
    #--------------------------------------------
    Hpol   = ꕕ(Had, Iₚ_1)                    # Matter
    Hpol  += ꕕ(Iₘ, Hₚ_1)                     # Photon
    Hpol  += ꕕ(µ, (â.T + â)) * χ_1           # Interaction
    Hpol  += ꕕ(µ @ µ, Iₚ_1) * (χ_1**2/(ωc/27.2114))    # Dipole Self-Energy

    #print ( "\tShape of Total Hamiltonian:", np.shape(Hpol) )

    print(f"\tMemory size of numpy array in (MB, GB): ({round(Hpol.size * Hpol.itemsize * 10 ** -6,2)},{round(Hpol.size * Hpol.itemsize * 10 ** -9,2)})" )

    return Hpol

def getAdiab(d,params):

    Nad = params.Nad # For all, set to large number (e.g., 500)

    E_POLARIZATION = np.array([ ( d == 'x' ), ( d == 'y' ), ( d == 'z' ) ])
    #print ("\tLight Polarization:", E_POLARIZATION * 1.0 )

    Had = np.zeros(( Nad + 1, Nad + 1 )) # Need to add one for ground state
    counter = 0
    for line in np.loadtxt(f"../TD-DFT/adiabtic_energies.dat").astype(float):
        state = int(line[0])
        energy = line[1]
        if ( state <= Nad ):
            Had[ counter, counter ] = energy / 27.2114 # Adiabatic Energies (convert to a.u.)
            #print ( state, energy )
            counter += 1

    dip_temp = np.loadtxt(f"../TD-DFT/dipole_matrix_E.dat") # 3-column file "i j (Dij)x (Dij)y (Dij)z". File length should be ~ Nad ** 2 from EOM-CCSD or TD-HF calcuation
    µ = np.zeros(( Nad + 1, Nad + 1 ))
    for line in dip_temp:
        i = int( line[0] )
        j = int( line[1] )
        if ( i > Nad or j > Nad ):
            continue
        µ[i,j] = np.dot( line[2:5], E_POLARIZATION ) # Choose 4 = mu_z, 3 = mu_y, 2 = mu_x, already in a.u.
        µ[j,i] = µ[i,j]

    #print ( np.shape(Had), np.shape(µ) )
    #print ( Had )
    #print ( µ )

    return Had, µ

def SolvePlotandSave(Hpol,Had,d,params):

        # Diagonalize polaritonic Hamiltonian and save
        E, U = SC_EXACT_EIGH(Hpol) # This is exact solution
        np.savetxt( f"data_PF/Epol_E{d}_eta{params.η_1}_wc{params.ωc}_Nad{params.Nad}.dat", E * 27.2114 )
        np.savetxt( f"data_PF/Upol_E{d}_eta{params.η_1}_wc{params.ωc}_Nad{params.Nad}.dat", U ) # These can be large
        #np.save( f"data_PF/Upol_E{d}_eta{params.η_1}_wc{params.ωc}_Nad{params.Nad}.dat", U ) # Smaller
        np.savetxt( f"data_PF/Char_E{d}_eta{params.η_1}_wc{params.ωc}_Nad{params.Nad}.dat", U[1,:] ** 2 ) # Photonic Character -- Saves Space Compared to Storing all U

        # Write original E_adiab to file for comparison
        np.savetxt( f"data_PF/Had.dat", Had[np.diag_indices(len(Had))] * 27.2114 )

def main_Serial():

    for d in ['z']: # CHOOSE POLARIZATION DIRECTION OF DIPOLE MOMENT
        Had, µ = getAdiab( d, params ) 
        Hpol = getĤpol( Had, µ, params )
        SolvePlotandSave( Hpol, Had, d, params)

    


if __name__ == "__main__":

    print ("\tParameters (η_1, wc):", params.η_1, params.ωc)
    main_Serial()