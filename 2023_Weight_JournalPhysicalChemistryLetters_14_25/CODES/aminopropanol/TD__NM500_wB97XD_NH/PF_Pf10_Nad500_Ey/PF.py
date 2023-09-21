import numpy as np
from scipy.linalg import eigh as SC_EXACT_EIGH
from numpy import kron as ꕕ
import sys
from matplotlib import pyplot as plt
import multiprocessing as mp
import subprocess as sp

def get_op_a(nf):
    a = np.zeros((nf,nf))
    for m in range(1,nf):
        a[m,m-1] = np.sqrt(m)
    return a.T

class params:

    Nad = 500 # Number of Electronic Basis States
    nf = 10 # Number of Fock Basis States
    
    η_1 = float( sys.argv[1] )
    ωc  = 3.0000 # eV

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
    â = get_op_a(nf) 
    #--------------------------------------------
    #       matter ⊗ photon 
    #--------------------------------------------
    Hpol   = ꕕ(Had, Iₚ_1)                    # Matter
    Hpol  += ꕕ(Iₘ, Hₚ_1)                     # Photon
    Hpol  += ꕕ(µ, (â.T + â)) * χ_1           # Interaction
    Hpol  += ꕕ(µ @ µ, Iₚ_1) * (χ_1**2/(ωc/27.2114))    # Dipole Self-Energy

    # BELOW IS ADDITION FOR TWO-MODE CAVITY
    #χ_3 = params.χ_3
    #Iₚ_3 = np.identity(nf)
    #Hₚ_3 = np.identity(nf)
    #Hₚ_3[np.diag_indices(nf)] = np.arange(nf) * (ωc/27.2114) * 3.0
    
    # Hpol   = ꕕ(ꕕ(Had, Iₚ_1),Iₚ_3)     # Matter
    # Hpol  += ꕕ(ꕕ(Iₘ, Hₚ_1),Iₚ_3)     # Photon 1
    # Hpol  += ꕕ(ꕕ(Iₘ, Iₚ_1),Hₚ_3)     # Photon 3
    # Hpol  += ꕕ(ꕕ(µ, (â.T + â)) * χ_1,Iₚ_3)      # Interaction 1
    # Hpol  += ꕕ(ꕕ(µ, Iₚ_1), (â.T + â)) * χ_3     # Interaction 3
    # Hpol  += ꕕ(ꕕ(µ @ µ, Iₚ_1),Iₚ_3) * (χ_1**2/(ωc/27.2114))    # Dipole Self-Energy
    # Hpol  += ꕕ(ꕕ(µ @ µ, Iₚ_1),Iₚ_3) * (χ_3**2/(3.0 * (ωc/27.2114))) # DSE 3

    #print ( "\tShape of Total Hamiltonian:", np.shape(Hpol) )


    print(f"\tMemory size of numpy array in (MB, GB): ({round(Hpol.size * Hpol.itemsize * 10 ** -6,2)},{round(Hpol.size * Hpol.itemsize * 10 ** -9,2)})" )

    return Hpol

def getAdiab(d,params):

    Nad = params.Nad # For all, set to large number (e.g., 500)

    E_POLARIZATION = np.array([ ( d == 'x' ), ( d == 'y' ), ( d == 'z' ) ])
    #print ("\tLight Polarization:", E_POLARIZATION * 1.0 )

    Had = np.zeros(( Nad + 1, Nad + 1 )) # Need to add one for ground state
    counter = 0
    for line in np.loadtxt(f"../adiabtic_energies.dat").astype(float):
        state = int(line[0])
        energy = line[1]
        if ( state <= Nad ):
            Had[ counter, counter ] = energy / 27.2114 # Adiabatic Energies (convert to a.u.)
            #print ( state, energy )
            counter += 1

    dip_temp = np.loadtxt(f"../dipole_matrix_E.dat") # 3-column file "i j (Dij)x (Dij)y (Dij)z". File length should be ~ Nad ** 2 from EOM-CCSD or TD-HF calcuation
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

def get_average_photon_number(U, NMax):
    photon_number = np.zeros(( len(U) ))
    I_m = np.identity( params.Nad + 1 )
    op_a = get_op_a( params.nf )
    N_Total = ꕕ( I_m, op_a.T @ op_a )
    for j in range( NMax ):
        photon_number[j] = np.conjugate(U[:,j]) @ N_Total @ U[:,j]
    return photon_number

def get_average_matter_number(U, NMax):
    matter_number = np.zeros(( len(U) ))
    I_ph = np.identity( params.nf )
    op_a = get_op_a( params.Nad+1 )
    N_Total = ꕕ( op_a.T @ op_a, I_ph )
    for j in range( NMax ):
        matter_number[j] = np.conjugate(U[:,j]) @ N_Total @ U[:,j]
    return matter_number



def SolvePlotandSave(Hpol,Had,d,params):

        # Diagonalize polaritonic Hamiltonian and save
        E, U = SC_EXACT_EIGH(Hpol) # This is exact solution
        np.savetxt( f"data_polariton/Epol_E{d}_eta{params.η_1}_wc{params.ωc}_Nad{params.Nad}.dat", E * 27.2114 )
        #np.savetxt( f"data_polariton/Upol_E{d}_eta{params.η_1}_wc{params.ωc}_Nad{params.Nad}.dat", U ) # These can be large
        #np.save( f"data_polariton/Upol_E{d}_eta{params.η_1}_wc{params.ωc}_Nad{params.Nad}.dat", U ) # Smaller
        #np.savetxt( f"data_polariton/Char_E{d}_eta{params.η_1}_wc{params.ωc}_Nad{params.Nad}.dat", U[1,:] ** 2 ) # Photonic Character -- Saves Space Compared to Storing all U
        np.savetxt( f"data_polariton/Photon_Number_E{d}_eta{params.η_1}_wc{params.ωc}_Nad{params.Nad}.dat", get_average_photon_number(U,50) ) # Photonic Number -- Saves Space Compared to Storing all U
        np.savetxt( f"data_polariton/Matter_Number_E{d}_eta{params.η_1}_wc{params.ωc}_Nad{params.Nad}.dat", get_average_matter_number(U,50) ) # Photonic Number -- Saves Space Compared to Storing all U


        # Write original E_adiab to file for comparison
        np.savetxt( f"data_polariton/Had.dat", Had[np.diag_indices(len(Had))] * 27.2114 )

def main_Serial():

    for d in ['y']: # CHOOSE POLARIZATION DIRECTION OF DIPOLE MOMENT
        Had, µ = getAdiab( d, params ) 
        Hpol = getĤpol( Had, µ, params )
        SolvePlotandSave( Hpol, Had, d, params)

   


if __name__ == "__main__":

    sp.call("mkdir data_polariton/",shell=True)

    print ("\tParameters (η_1, wc):", params.η_1, params.ωc)
    main_Serial()