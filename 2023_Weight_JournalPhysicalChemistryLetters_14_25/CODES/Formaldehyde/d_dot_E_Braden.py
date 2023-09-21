import numpy as np
from scipy.linalg import eigh as SC_EXACT_EIGH
from scipy.sparse import linalg as SC_LANCZOS
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

    startFolder = int(sys.argv[1]) # Folder name which houses the adiabatic electronic structure energies and dipole moments [get_HAM_and_DIP_Matrix.py]
    endFolder = int(sys.argv[2])

    NCPUS = 24 # Parallelize over {NCPUS} cores
    Nad = 10 # Number of Electronic Basis States
    nf = 30 # Number of Fock Basis States
    
    # Flick and Narang: Panel (b)
    #η_1  = 0.04 # Coupling Strength
    #ωc = 8.71/27.2114
    
    # Flick and Narang: Panel (c)
    #η_1  = 0.04 # Coupling Strength
    #ωc = 8.16/27.2114

    η_1 = float( sys.argv[3] )
    ωc  = float( sys.argv[4] ) # in eV

    χ_1 = (ωc/27.2114) * η_1
    
    # BELOW IS ADDITION FOR TWO-MODE CAVITY
    #η_3  = 0.15 / np.sqrt(3)
    #χ_3 = ωc * η_3 * 3

def getĤpol(Had, µ, params):
    """
    Input: Had, adiabatic hamiltonian (diagonal) energies from electronic structure
    Output: Hpol, Hamiltonian after Kron. product with single- or double-mode photon cavity space
    """

    ns = len(Had)
    nf = params.nf
    ωc = params.ωc 
    χ_1 = params.χ_1
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

    return Hpol

def getAdiab(d,folderNAME,params):

    Nad = params.Nad # For all, set to large number (e.g., 500)

    E_POLARIZATION = np.array([ ( d == 'x' ), ( d == 'y' ), ( d == 'z' ) ])
    #print ("\tLight Polarization:", E_POLARIZATION * 1.0 )

    SYMM = [ j.split("-")[1].split("\n")[0] for j in open(f"{folderNAME}/all_Exc_Sym.dat", "r").readlines() ] # Line = "Singlet-A1\n"
    A1States = [ j+1 for j in range(len(SYMM)) if SYMM[j] == "A1" and j < Nad  ] # These are excited states, so the index must be shifted such that S1 = 1
    A1States = np.append( np.array([0]), np.array(A1States) ) # Add ground state
    NA1States = len(A1States)
    #print ("\tStates with A1 Symmetry: S_k =", A1States)
    #print (NA1States)

    Had = np.zeros(( NA1States, NA1States )) # Need to add one for ground state
    counter = 0
    for line in np.loadtxt(f"{folderNAME}/adiabtic_energies.dat").astype(float):
        state = int(line[0])
        energy = line[1]
        if ( state < Nad and state in A1States ):
            Had[ counter, counter ] = energy / 27.2114 # Adiabatic Energies (convert to a.u.)
            #print ( state, energy )
            counter += 1

    dip_temp = np.loadtxt(f"{folderNAME}/dipole_matrix_E.dat") # 3-column file "i j (Dij)x (Dij)y (Dij)z". File length should be ~ Nad ** 2 from EOM-CCSD or TD-HF calcuation
    µ = np.zeros(( NA1States, NA1States ))
    for line in dip_temp:
        i = int( line[0] )
        j = int( line[1] )
        if ( i >= Nad or j >= Nad or i not in A1States or j not in A1States ): # Skip reading dipoles if we are not including them.
            continue
        indi = np.where( A1States == i )[0][0]
        indj = np.where( A1States == j )[0][0]
        #print (i, j, np.where( A1States == i )[0][0], np.where( A1States == j )[0][0], np.dot( line[2:5], E_POLARIZATION ) )
        µ[indi,indj] = np.dot( line[2:5], E_POLARIZATION ) # Choose 4 = mu_z, 3 = mu_y, 2 = mu_x, already in a.u.
        µ[indj,indi] = µ[indi,indj]

    #print ( np.shape(Had), np.shape(µ) )
    #print ( Had )
    #print ( µ )

    return Had, µ

def SolvePlotandSave(Hpol,Had,d,folderNAME,params):

        def get_photon_number( U, a, Had ):
            Nm = len(Had)
            Nf = 30
            photon_number = np.zeros(( len(U) ))
            I_m = np.identity( Nm )
            adaga_total = ꕕ( I_m, a.T @ a )
            #print( np.shape(adaga_total), len(U) )
            for j in range( len(U) ):
                photon_number[j] = U[:,j].T @ adaga_total @ U[:,j]
            return photon_number

        # Diagonalize polaritonic Hamiltonian and save
        #E,U = np.real( SC_LANCZOS.eigs(Hpol,k=10) ) # This uses the Lanczos routines to find lowest k modes
        E, U = SC_EXACT_EIGH(Hpol) # This is exact solution
        np.savetxt( f"pol_data/Epol_E{d}_A0{np.round(params.η_1,4)}_wc{np.round(params.ωc,4)}_{folderNAME}.dat", E * 27.2114 )
        #np.savetxt( f"Upol_E{d}_{folderNAME}.dat", U ) # These can be large
        #photon_number = get_photon_number( U, ĉ(30), Had )
        #np.savetxt( f"pol_data/photon_number_E{d}_A0{np.round(params.η_1,4)}_wc{np.round(params.ωc,4)}_{folderNAME}.dat", photon_number )

        # Plot Hamiltonian elements
        """
        plt.contourf( np.log( np.abs(Hpol) + 0.01 ) )
        plt.colorbar()
        plt.title( "Log [Ĥ_Pol. (a.u.)]" )
        plt.savefig(f"Hpol_E{d}_{folderNAME}.jpg")
        plt.clf()
        """
        # Write original E_adiab to file for comparison
        np.savetxt( f"pol_data/Had_E{d}_{folderNAME}.dat", Had[np.diag_indices(len(Had))] * 27.2114 )

def main_Serial():

    for d in ['z']: # CHOOSE POLARIZATION DIRECTION OF DIPOLE MOMENT
        for folderNAME in range( params.startFolder, params.endFolder+1 ):
            print (f"Step: {folderNAME}")

            Had, µ = getAdiab( d, folderNAME, params ) 
            #print ( f"\tShape of Total Hamiltonian: ({ len(Had) * params.nf } x { len(Had) * params.nf }) " )
            Hpol = getĤpol( Had, µ, params )
            SolvePlotandSave( Hpol, Had, d, folderNAME, params)

def main_Parallel(folderNAME):

    for d in ['z']: # CHOOSE POLARIZATION DIRECTION OF DIPOLE MOMENT
        
        #print (f"Step: {folderNAME}")
        Had, µ = getAdiab( d, folderNAME, params ) 
        #print ( f"\tShape of Total Hamiltonian: ({ len(Had) * params.nf } x { len(Had) * params.nf }) " )
        Hpol = getĤpol( Had, µ, params )
        SolvePlotandSave( Hpol, Had, d, folderNAME, params)

    


if __name__ == "__main__":

    print ("\tParameters (η_1, wc):", params.η_1, params.ωc)

    runList = np.arange( params.startFolder, params.endFolder+1 )
    with mp.Pool(processes=params.NCPUS) as pool:
        pool.map(main_Parallel,runList)
    #main_Serial()