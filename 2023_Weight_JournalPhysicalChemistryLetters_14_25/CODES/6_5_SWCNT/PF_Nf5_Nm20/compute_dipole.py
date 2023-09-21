import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from numba import jit
#from pygifsicle import optimize as gifOPT # This needs to be installed somewhere
#from PIL import Image, ImageDraw, ImageFont
#import imageio
import subprocess as sp
from time import time, sleep
from scipy.integrate import simps
from scipy.interpolate import interp2d
from scipy.ndimage import gaussian_filter
from scipy.special import hermite
from scipy.special import factorial
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy.interpolate import interp1d



"""
Install pygifcicle:
pip3 install pygifsicle


Install gifsicle: ( echo "$(pwd)" = /scratch/bweight/software/ )
curl -sL http://www.lcdf.org/gifsicle/gifsicle-1.91.tar.gz | tar -zx
cd gifsicle-1.91
./configure --disable-gifview
make install exec_prefix=$(pwd) prefix=$(pwd) datarootdir=$(pwd)
"""

def getGlobals():
    global eta_list, wc, Nad, d, NPolCompute, Neta, write_TD_Files, data_PF, print_level
    global CMomFile, DMomFile, QMomFile, QGrid, plotDIM
    #eta_list = [ 0.0, 0.01, 0.02 ]  # np.array([ 0.0, 0.01, 0.05, 0.1, 0.2 ])  #np.arange( 0.0, 0.2+0.05, 0.05 )
    eta_list = np.arange( 0.0, 0.02+0.0005, 0.0005 )
    Neta = len(eta_list)
    wc = 2.0 # 1.974 # 1.9733
    Nad = 20
    NPolCompute = 30
    QGrid = np.arange( -20, 20, 0.2 )
    write_TD_Files = True
    print_level = 0 # 0 = Minimal, 1 = Some, 2 = Debugging
    d = 'z'
    plotDIM = 'z' # Dimension for plotting 2D transition density
    data_PF = "data_PF"
    sp.call('mkdir -p data_expansion',shell=True)
    sp.call('mkdir -p data_dipole',shell=True)
    sp.call('mkdir -p data_TD',shell=True)
    sp.call('mkdir -p data_ElectricMoments',shell=True)
    sp.call('mkdir -p data_spectra',shell=True)


    # Instantiate Files
    CMomFile = open("data_ElectricMoments/Exciton_Transition_Charge.dat","w")
    DMomFile = open("data_ElectricMoments/Exciton_Transition_Dipole.dat","w")
    QMomFile = open("data_ElectricMoments/Exciton_Transition_Quadrupole.dat","w")
    CMomFile.write(f"# j k Charge\n")
    DMomFile.write(f"# j k Dx Dy Dz\n")
    QMomFile.write(f"# j k Qxx Qyy Qzz Qxy Qxz Qyz\n")



#@jit(nopython=True)
def get_TD_fast_1r( Upol, TD_matter, NPolCompute ):
    # EXACT:  <P0| TD(R) |Pa> = SUM_{j k n} C_{jn}^(0) C_{kn}^(a) TD(R)_{jk}
    # Note: We only have T_{00} and T_{0k} matrix elements
    TD_pol = np.zeros(( len(Upol), NPolCompute, Nxyz[0]*Nxyz[1]*Nxyz[2] )) # Eta, NPol, NGrid, assuming calculation of P0 to Pa, a = 0,1,2
    for etaIND in range(Neta):
        for a in range( NPolCompute ):                       # Choose only up to what is needed..
            for n in range( NFock ):                         # Photons
                for state_j in range( Nad+1 ):               # Matter
                    polIND1 = state_j * NFock  + n           # Only all | 0,n > basis here
                    for state_k in range( Nad+1 ):           # Matter
                        polIND2 = state_k * NFock  + n       # Including all | S,n > basis here
                        TD_pol[etaIND,a,:] += Upol[etaIND,polIND1,0] * Upol[etaIND,polIND2,a] * TD_matter[state_j,state_k,:,:,:].flatten()
    return TD_pol

def getHO_q( n ):
    print ( f"HO_q_{n}" ,n )
    wc_au = wc / 27.2114
    HO_q_n = 1/np.sqrt(2**n * factorial(n) ) * ( wc_au/np.pi )**(1/4) * np.exp( - wc_au * QGrid**2 / 2 ) * hermite(n)( np.sqrt(wc_au) * QGrid)
    HO_q_n /= np.max( HO_q_n )

    plt.plot( QGrid, HO_q_n )
    plt.xlabel( r'q$_c$ (a.u.)', fontsize=15 )
    plt.ylabel( r'$\phi (q_c)$', fontsize=15 )
    plt.tight_layout()
    plt.savefig(f"HO_q_{n}.jpg" )
    plt.clf()

    return HO_q_n 

def getTraceTD( TD_matter ):
    dict_xyz = { 'x':(1+2,2+2), 'y':(0+2,2+2), 'z':(0+2,1+2) } # Adding two for state_j and state_k labels
    TD_projected = np.sum( TD_matter, axis=dict_xyz[plotDIM] )
    return TD_projected

@jit(nopython=True)
def get_TD_fast_1r1q( Upol, TD_projected, NPolCompute, HO_WFNs  ):
    # EXACT:  <P0| TD(R,q) |Pa> = SUM_{j k n m} C_{jn}^(0) C_{km}^(a) TD(R)_{jk} \phi(q)_n^HO \phi(q)_m^HO
    # Note: We only have T_{00} and T_{0k} matrix elements

    TD_pol = np.zeros(( len(Upol), NPolCompute, np.shape(TD_projected)[-1],  len(QGrid) )) # Eta, NPol, RGrid, QGrid ... assuming calculation of P0 to Pa, a = 0,1,2
    for etaIND in range(Neta):
        #for a in range( NPolCompute ):                       # Choose only up to what is needed..
        for a in [ 0, 1, 6, 7, 21, 22 ]:                                          # Choose only up to what is needed..
            print ( a )
            for n in range( NFock ):                         # Photon
                for state_j in range( Nad+1 ):               # Matter
                    polIND1 = state_j * NFock  + n               # Only all | S_j,n > basis here
                    for state_k in range( Nad+1 ):           # Matter
                        for m in range( NFock ):             # Photon
                            polIND2 = state_k * NFock  + m       # Including all | S_k,n > basis here
                            for qi in range( len(QGrid) ):
                                TD_pol[etaIND,a,:,qi] += Upol[etaIND,polIND1,0] * Upol[etaIND,polIND2,a] * TD_projected[state_j,state_k,:] * (HO_WFNs[n])[qi] * (HO_WFNs[m])[qi]
    return TD_pol #.reshape(( len(Upol), NPolCompute, np.shape(TD_projected)[-1]  , len(QGrid) ))


#@jit(nopython=True)
def compute_dipole_1_fast( NFock, Upol, MU, NPolCompute ):
    # <P0| µ |Pa> = SUM_{j j' n} C_{in}^(0) C_{jn}^(a) µ_{ij}
    eff_dipole_1 = np.zeros(( len(Upol), NPolCompute ))
    for a in range( NPolCompute ): # This will give ground to excited state "a" dipole matrix elements. Choose only up to what is needed...
        for n in range( NFock ):                         # Photon
            for i in range( len(Upol[0])//NFock ):                   # Matter
                polIND1 = i *NFock  + n
                for j in range( len(Upol[0])//NFock ):               # Matter
                    polIND2 = j *NFock  + n
                    eff_dipole_1[:,a] += Upol[:,polIND1,0] * Upol[:,polIND2,a] * MU[i,j] # Implicit loop over eta
    return eff_dipole_1

#@jit(nopython=True)
def compute_dipole_2_fast( NFock, Upol, MU, NFockInclude, NMatterInclude, NPolCompute ):
    # |<P0| µ |Pa>|^2 = SUM_{j j' n} SUM_{ k k', m } C_{in}^(0) C_{i'm}^(0) C_{jn}^(a) C_{j'm}^(a) µ_{ij} µ_{i'j'}  (Someone should check this ~ BMW)
    eff_dipole_2 = np.zeros(( len(Upol), NPolCompute ))
    for a in range( NPolCompute ): # This will give ground to excited state "a" dipole matrix elements squared. Choose only up to what is needed...
        for n in range( NFockInclude ):                         # Photon
            for i in range( NMatterInclude ):                   # Matter
                for ip in range( NMatterInclude ):              # Matter
                    for m in range( NFockInclude ):             # Photon
                        for j in range( NMatterInclude ):       # Matter
                            for jp in range( NMatterInclude ):  # Matter
                                polIND1  = i *NFock  + n            
                                polIND1p = ip*NFock  + m
                                polIND2  = j *NFock  + n
                                polIND2p = jp*NFock  + m
                                
                                eff_dipole_2[:,a] += Upol[:,polIND1,0] * Upol[:,polIND1p,0] * Upol[:,polIND2,a] * Upol[:,polIND2p,a] * MU[i,j] * MU[ip,jp] # Implicit loop over eta
    return eff_dipole_2

#@jit(nopython=True)
def compute_eff_osc_str_fast( Npolar, Epol, eff_dipole_1 ):
    eff_osc_str = np.zeros(( len(eff_dipole_1), Npolar ))
    for a in range( NPolCompute ):
        dE = (Epol[:,a] - Epol[:,0]) / 27.2114 # eV to a.u.
        eff_osc_str[:,a] = (2/3) * dE[:] * eff_dipole_1[:,a] ** 2 # Implicit loop over eta
    return eff_osc_str

#@jit(nopython=True)
def trace_fock( Upol ):
    Upol_ad = np.zeros(( len(eta_list), NPolar, Nad + 1 ))      # Coupling Strength, # of Polaritons, # of Adiabatic states
    for j in range( Nad + 1):                               # LOOP OVER MATTER BASIS
        for n in range( NFock ):                            # LOOP OVER FOCK BASIS
            polIND1 = j * NFock  + n                            # POL. SUPER-INDEX 1, DEFINES STATE OF INTEREST
            for jp in range( Nad + 1 ):                     # LOOP OVER MATTER BASIS
                for npp in range( NFock ):                  # LOOP OVER FOCK BASIS
                    polIND2 = jp * NFock + npp                  # POL. SUPER-INDEX 2, DEFINES EXPANSION COMPONENT
                    Upol_ad[ :, polIND1, jp ] += Upol[ :, polIND2, polIND1 ] ** 2 # (g,0), (g,1), ..., (g,NFock-1), (e,0), (e,1), ..., 
    return Upol_ad

#@jit(nopython=True)
def trace_ad( Upol ):
    Upol_fock = np.zeros(( len(eta_list), NPolar, NFock ))
    for j in range( Nad + 1 ):
        for n in range( NFock ):
            polIND1 = j * NFock  + n
            for jp in range( Nad + 1 ):
                for npp in range( NFock ):
                    polIND2 = jp * NFock + npp
                    Upol_fock[ :, polIND1, npp ] += Upol[ :, polIND2, polIND1 ] ** 2 # (g,0), (g,1), ..., (e,0), (e,1), ..., 
    return Upol_fock

def get_HadMU():
    global NFock, NPolar

    # Get Data
    NPolar = len(np.loadtxt(f"{data_PF}/Epol_E{d}_eta{eta_list[0]}_wc{wc}_Nad{Nad}.dat"))
    Epol = np.zeros(( len(eta_list), NPolar ))
    Char = np.zeros(( len(eta_list), NPolar ))
    Upol = np.zeros(( len(eta_list), NPolar, NPolar ))
    for j in range( len(eta_list) ):
        eta = round( eta_list[j], 8 )
        print ( f"Reading Upol file {j} of {len(eta_list)}" )
        Epol[j,:] = np.loadtxt(f"{data_PF}/Epol_E{d}_eta{eta}_wc{wc}_Nad{Nad}.dat")
        Char[j,:] = np.loadtxt(f"{data_PF}/Char_E{d}_eta{eta}_wc{wc}_Nad{Nad}.dat")
        #Upol[j,:,:] = np.load(f"data/Upol_E{d}_eta{eta}_wc{wc}_Nad{Nad}.dat.npy")
        Upol[j,:,:] = np.loadtxt(f"{data_PF}/Upol_E{d}_eta{eta}_wc{wc}_Nad{Nad}.dat")
     
    NFock = len( Upol[0,0,:] ) // (Nad+1)
    print (f"NFock = {NFock}")

    # Get molecular dipoles
    E_POLARIZATION = np.array([ ( d == 'x' ), ( d == 'y' ), ( d == 'z' ) ])
    dip_temp = np.loadtxt(f"../TD-DFT/dipole_matrix_E.dat") # 3-column file "i j (Dij)x (Dij)y (Dij)z". File length should be ~ Nad ** 2 from EOM-CCSD or TD-HF calcuation
    MU = np.zeros(( Nad + 1, Nad + 1 ))
    for line in dip_temp:
        i = int( line[0] )
        j = int( line[1] )
        if ( i > Nad or j > Nad ):
            continue
        MU[i,j] = np.dot( line[2:5], E_POLARIZATION )
        MU[j,i] = MU[i,j]
    
    assert( NPolCompute <= len(Upol[0]) ), "Number of requested polaritons to compute is larger than the number of polaritons."

    return Epol, Char, Upol, MU

def compute_Electrostatic_Moments( TD, state_j, state_k ):
    """
    Input:
        TD: nd.array(( NStates,NStates,Nx,Ny,Nz )) - transition density (position-basis)
        state_j: int - electronic state j
        state_k: int - electronic state k
    Returns:
        None
    """

    xyz_dic = { (0,1):2, (0,2):1, (1,2):0 }

    # Charge
    charge = np.sum(TD[state_j,state_k,:,:,:])*dLxyz[0]*dLxyz[1]*dLxyz[2]
    CMomFile.write( f"{state_j}\t{state_k}\t{charge}\n" )
    if ( print_level >= 1 ):
        print (f'\t\t,CHARGE of T_{state_j}-{state_k} = {charge} =? 0.000')

    # Dipole
    op_Rx = np.arange(0,Nxyz[0])*dLxyz[0]
    op_Ry = np.arange(0,Nxyz[1])*dLxyz[1]
    op_Rz = np.arange(0,Nxyz[2])*dLxyz[2]
    op_R = [ op_Rx, op_Ry, op_Rz ]
    

    Dx = np.sum( np.sum(TD[state_j,state_k,:,:,:],axis=(1,2)) * op_R[0] )*dLxyz[0]*dLxyz[1]*dLxyz[2]
    Dy = np.sum( np.sum(TD[state_j,state_k,:,:,:],axis=(0,2)) * op_R[1] )*dLxyz[0]*dLxyz[1]*dLxyz[2]
    Dz = np.sum( np.sum(TD[state_j,state_k,:,:,:],axis=(0,1)) * op_R[2] )*dLxyz[0]*dLxyz[1]*dLxyz[2]

    DMomFile.write( f"{state_j}\t{state_k}\t{Dx}\t{Dy}\t{Dz}\n" )

    if ( print_level >= 1 ):
        print (f'\t\t\tDipole Moment (X) of T_{state_j}-{state_k} = {Dx}')
        print (f'\t\t\tDipole Moment (Y) of T_{state_j}-{state_k} = {Dy}')
        print (f'\t\t\tDipole Moment (Z) of T_{state_j}-{state_k} = {Dz}')
    
    # Quadrupole, depends on location of origin. Shift coordinates to Q_COM
    # V_QUAD = 1/(4 \pi \eps_o) 1/(2 r^3) SUM_[jk] \hat{r}\hat{j} Q_{jk}
    # Q_{jk} = INT ( 3 r_j r_k - r_j^2 (j==k) ) \rho(r)  dr d^3r
    Q_COM = [ np.average(op_R[0]), np.average(op_R[1]), np.average(op_R[2])  ] # Find center of box
    op_R[0] -= Q_COM[0] # Shift coordinates to center of box
    op_R[1] -= Q_COM[1] # Shift coordinates to center of box
    op_R[2] -= Q_COM[2] # Shift coordinates to center of box
    op_Q = np.zeros(( 3, 3 ))
    for d1 in range( 3 ): # x,y,z
        for d2 in range( d1, 3 ): # x,y,z
            
            if ( d1 == d2 ):
                ind_1 = (d1+1)%3
                ind_2 = (d1+2)%3
                T_proj = np.sum(TD[state_j,state_k,:,:,:],axis=(ind_1, ind_2)) * dLxyz[ind_1] * dLxyz[ind_2]
                #op_Q[d1,d2] = np.sum( np.sum( T_proj[:] * op_R[d1]**2,axis=0 ) * op_R[d2] )*dLxyz[d1]
                op_Q[d1,d2] = 2 * np.sum( T_proj[:] * op_R[d1]**2 ) * dLxyz[d1]

            else:
                dTr = xyz_dic[ (d1,d2) ]
                T_proj = np.sum( TD[state_j,state_k,:,:,:],axis=(dTr) ) * dLxyz[dTr]
                tmp = []
                for j in range( Nxyz[d2] ):
                    tmp.append( 3 * np.sum( T_proj[:,d2] * op_R[d1]**2, axis=0 ) * dLxyz[d1] )
                op_Q[d1,d2] = np.sum( np.array(tmp) * op_R[d2] )*dLxyz[d2]
                op_Q[d2,d1] = op_Q[d1,d2]

    # j k Qxx Qyy Qzz Qxy Qxz Qyz
    QMomFile.write( f"{state_j}\t{state_k}\t{op_Q[0,0]}\t{op_Q[1,1]}\t{op_Q[2,2]}\t{op_Q[0,1]}\t{op_Q[0,2]}\t{op_Q[1,2]}\n" )

    if ( print_level >= 1 ):
        print ( "Octopole Tensor\n", op_Q )
        print ( 'Tr[Q] =', np.round( np.sum(op_Q[np.diag_indices(3)]) ,4) )

def get_TD_Data():
    """
    NOTE: WE ARE ONLY CALCULATING S0 to Sk TRANSITION DENSITY. IN PRINCIPLE, WE NEED ALL TERMS.
    """

    #print ("\tStarting to Read TD Files.")

    # Get size from first TD file
    global NAtoms, NGrid, Nxyz, dLxyz, Lxyz, coords
    header = np.array([ j.split() for j in open(f"../TD-DFT/TDMat/trans-0_1.cube","r").readlines()[1:6] ])
    NAtoms = int(header[1,0])
    NGrid  = int( header[0,1] )
    Nxyz   = (header[2:,0]).astype(int)
    dLxyz  = np.array([header[2,1],header[3,2],header[4,3] ]).astype(float)
    #Lxyz   = np.array([ Nxyz[0]*dLxyz[0], Nxyz[1]*dLxyz[1], Nxyz[2]*dLxyz[2] ])
    Lxyz   = np.array([ header[1,1], header[1,2], header[1,3] ]).astype(float)
    if ( Lxyz[0] < 0 ): Lxyz *= -1.000 # Switch Sign, Already Angstroms
    if ( Lxyz[0] > 0 ): Lxyz *= 0.529 # Convert from Bohr to Angstroms
    Vol    = Lxyz[0] * Lxyz[1] * Lxyz[2]
    

    print (f'\tNAtoms   = {NAtoms}')
    print (f'\tTD Grid  = {NGrid}')
    print (f'\tNx Ny Nz = {Nxyz[0]} {Nxyz[1]} {Nxyz[2]}')
    print (f'\tLx Ly Lz = {Lxyz[0]} {Lxyz[1]} {Lxyz[2]} A')
    print (f'\tVolume   = {Vol} A^3')

    NStart = NAtoms + 6

    coords = np.array([ j for j in open(f"../TD-DFT/TDMat/trans-0_1.cube","r").readlines()[6:NStart] ])


    TD = np.zeros(( Nad+1, Nad+1, Nxyz[0], Nxyz[1], Nxyz[2]  ))
    print(f"\tMemory size of transition density array in (MB, GB): ({round(TD.size * TD.itemsize * 10 ** -6,2)},{round(TD.size * TD.itemsize * 10 ** -9,2)})" )
    for state_j in range(Nad+1):
        for state_k in range(state_j,Nad+1):
            #print (f'\tReading Transition Density {mat+1} of {Nad}.')
            temp = []
            try:
                lines = open(f"../TD-DFT/TDMat/trans-{state_j}_{state_k}.cube",'r').readlines()[NStart:]
            except FileNotFoundError:
                print (f'\t****** File "trans-{state_j}_{state_k}.cube" not found. Skipping this matrix element. ******')
                continue
            for count, line in enumerate(lines):
                t = line.split('\n')[0].split()
                for j in range(len(t)):
                    temp.append( float(t[j]) )
            TD[state_j,state_k,:,:,:] = np.array( temp ).reshape(( Nxyz[0],Nxyz[1],Nxyz[2] ))
            if ( state_j != state_k ):
                TD[state_k,state_j,:,:,:] = 1.00 * TD[state_j,state_k,:,:,:] # SHOULD THIS SYMMETRIC, RIGHT ???
            elif( state_j == state_k and state_j == 0 ):
                norm = np.sum( TD[state_j,state_k,:,:,:],axis=(0,1,2) ) * dLxyz[0]*dLxyz[1]*dLxyz[2]
                TD[state_j,state_k,:,:,:] /= norm # Normalize to 1.0 instead of NELECT
            compute_Electrostatic_Moments( TD, state_j, state_k )                   
                    

    for dim in ['x','y','z']:
        if ( dim == 'x' ):
            for state_j in [0,6]: #Nad+1):
                for state_k in np.arange(state_j,Nad+1):
                    if ( state_j == state_k and state_j != 0 ): continue
                    f = np.sum( TD[state_j,state_k,:,:,:],axis=(1,2) )
                    grid = np.arange(Nxyz[0]) * 0.529 * dLxyz[0]
                    plt.plot( grid, f, '-' , label=f'TD_{state_j}-{state_k}' )
        elif( dim == 'y' ):
            for state_j in [0,6]: #Nad+1):
                for state_k in np.arange(state_j,Nad+1):
                    if ( state_j == state_k and state_j != 0 ): continue
                    f = np.sum( TD[state_j,state_k,:,:,:],axis=(0,2) )
                    grid = np.arange(Nxyz[1]) * 0.529 * dLxyz[1]
                    plt.plot( grid, f, '-' , label=f'TD_{state_j}-{state_k}' )
        elif( dim == 'z' ):    
            for state_j in [0,6]: #Nad+1):
                for state_k in np.arange(state_j,Nad+1):
                    if ( state_j == state_k and state_j != 0 ): continue
                    f = np.sum( TD[state_j,state_k,:,:,:],axis=(0,1) )
                    grid = np.arange(Nxyz[2]) * 0.529 * dLxyz[2]
                    plt.plot( grid, f, '-' , label=f'TD_{state_j}-{state_k}' )


        plt.legend()
        plt.tight_layout()
        plt.savefig(f"data_TD/Transition_Density_Slices_{dim}.jpg")
        plt.clf()


    return TD

def getExansion(Upol):

    print ("Getting Projected Expansion Coefficients")

    Upol_fock = trace_ad(Upol) # Eta, Polariton, Adiabatic State
    Upol_ad = trace_fock(Upol) # Eta, Polariton, Fock State

    #print ("Check Total Traces:")
    #print (f"Upol_fock: { np.sum(Upol_fock[:,5,:],axis=1) }")
    #print (f"Upol_ad: { np.sum(Upol_ad[:,5,:],axis=1) }")
    assert( np.allclose(np.sum(Upol_fock[:,5,:],axis=1), np.ones(len(eta_list))) ), "Fock traces are not normalized to unity."
    assert( np.allclose(np.sum(Upol_ad[:,5,:],axis=1), np.ones(len(eta_list))) ), "Matter traces are not normalized to unity."

    # Plot Projections
    for etaIND in range(len(eta_list)):
        pLIST = [0,1,2,3,4,5,6,7,8,9,10]
        for p in pLIST:
            plt.plot( np.arange( NFock ), Upol_fock[etaIND, p, : ], label=f"P{p}" )
        plt.legend(loc='upper right')
        plt.xlabel('Fock State |n>',fontsize=12)
        plt.ylabel('Population',fontsize=12)
        plt.savefig(f'data_expansion/Upol_Fock_eta{round(eta_list[etaIND],8)}.jpg')
        plt.clf()

        pLIST = [0,1,2,3,4,5,6,7,8,9,10]
        for p in pLIST:
            plt.plot( np.arange( Nad + 1 ), Upol_ad[etaIND, p, : ], label=f"P{p}" )
        plt.legend(loc='upper right')
        plt.xlabel('Matter State |j>',fontsize=12)
        plt.ylabel('Population',fontsize=12)
        plt.savefig(f'data_expansion/Upol_Matter_eta{round(eta_list[etaIND],8)}.jpg')
        plt.clf()

        pLIST = [0,1,2,3,4,5,6,7,8,9,10]
        for p in pLIST:
            plt.plot( np.arange( NPolar ), Upol[etaIND, p, : ]**2, label=f"P{p}" )
        plt.legend(loc='upper right')
        plt.xlim(-25,1500)
        plt.xlabel('Basis State |j>|n>',fontsize=12)
        plt.ylabel('Population',fontsize=12)
        plt.savefig(f'data_expansion/Upol_Coeff_eta{round(eta_list[etaIND],8)}.jpg')
        plt.clf()

    # Plot population of Pj as function of eta
    pLIST = [0,1,2,3,4,5,6,7,8,9,10]
    ad_list = [0,1,2,3,4,5,6]
    for p in pLIST:
        for ad in ad_list:
            if( np.max(Upol_ad[:, p, ad ]**2) > 0.01 ):
                plt.plot( eta_list, Upol_ad[:, p, ad ]**2,"o-", label=f"S$_{ad}$ in P$_{p}$" )
        plt.xlabel('Coupling Strength, $A_0$',fontsize=12)
        plt.ylabel('Population',fontsize=12)
        plt.legend()
        plt.savefig(f'data_expansion/Upol_Matter_Pop_in_P{p}__wc{round(wc,3)}.jpg')
        plt.clf()


def getDipole_pow_1(Upol,MU):

    print("Starting ground to excited transition dipole moments.")

    eff_dipole_1 = np.zeros(( len(eta_list), NPolCompute))

    eff_dipole_1[:,:] = compute_dipole_1_fast( NFock, Upol, MU, NPolCompute )
    
    """
    # Plot Dipole Convergence
    outside_cavity = [ 0.0000,  0.0000, 1.23870, 0.00020, 1.87960, 0.00260, 26.35430, 0.00160 ] # 00, 01, 02, 03, 04, 05
    
    for etaIND in range(len(eta_list)):
        
        pLIST = [0,1,2,3,4,5,6,7,8,9,10]
        for p in pLIST:
            plt.plot( np.arange( (Nad+1)//1 ), np.abs( eff_dipole_1[etaIND, p, -1, :(Nad+1)//1 ] ), label=f"P{p}" ) # Choose Most Fock States
        for ex in range(8):
            plt.plot( np.arange( Nad+1 ), outside_cavity[ex]*np.ones( Nad+1 ), '--', c='black'  )
        plt.legend(loc='upper right')
        plt.xlim(0,Nad)
        plt.xlabel('Matter State |j>',fontsize=12)
        plt.ylabel('Ground-to-Excited Transition Dipole (0-->a)',fontsize=12)
        plt.tight_layout()
        plt.savefig(f'data_dipole/Dipole_MatterSCAN_eta{round(eta_list[etaIND],8)}.jpg')
        plt.clf()

        pLIST = [0,1,2,3,4,5,6,7,8,9,10]
        for p in pLIST:
            plt.plot( np.arange( NFock ), np.abs( eff_dipole_1[etaIND, p, :, (Nad+1)//1-1 ] ), label=f"P{p}" ) # Choose Most Matter States
        for ex in range(8):
            plt.plot( np.arange( NFock ), outside_cavity[ex]*np.ones( NFock ), '--', c='black'  )
        plt.legend(loc='upper right')
        plt.xlim(0,NFock-1)
        plt.xlabel('Fock State |n>',fontsize=12)
        plt.ylabel('Ground-to-Excited Transition Dipole (0-->a)',fontsize=12)
        plt.tight_layout()
        plt.savefig(f'data_dipole/Dipole_FockSCAN_eta{round(eta_list[etaIND],8)}.jpg')
        plt.clf()


        pLIST = [ j for j in range(NPolCompute) ]
        for p in pLIST:
            if ( np.max(eff_dipole_1[:, p, -1, -1 ]) > 3 ):
                plt.plot( eta_list, np.abs( eff_dipole_1[:, p, -1, -1 ] ), label=f"P{p}" ) # Choose Most Matter/Fock States
            else:
                plt.plot( eta_list, np.abs( eff_dipole_1[:, p, -1, -1 ] ) ) # Choose Most Matter/Fock States
        for ex in range(8):
            plt.plot( eta_list, outside_cavity[ex]*np.ones( len(eta_list) ), '--', c='black'  )
        plt.legend(loc='upper right')
        plt.xlim(0,eta_list[-1])
        plt.xlabel('Coupling Strength A$_o$',fontsize=12)
        plt.ylabel('Ground-to-Excited Transition Dipole (0-->a)',fontsize=12)
        plt.tight_layout()
        plt.savefig(f'data_dipole/Dipole_Polaritons_Coupling.jpg')
        plt.clf()
    """
    # Print Results
    #for etaIND in range(len(eta_list)):
    #    for p in range( NPolCompute ):
    #        np.savetxt(f'data_dipole/dipole_pow_1__MatterSCAN_NFock{NFock}_eta{round(eta_list[etaIND],8)}_P{p}.dat', np.c_[ np.arange(Nad+1), eff_dipole_1[ etaIND, p, -1, : ] ] )
    #        np.savetxt(f'data_dipole/dipole_pow_1__FockSCAN_NMatter{Nad}_eta{round(eta_list[etaIND],8)}_P{p}.dat', np.c_[ np.arange(NFock), eff_dipole_1[ etaIND, p, :, -1 ] ] )

    for etaIND in range(len(eta_list)):
        dipFile = open( f'data_dipole/dipole_Polaritons_eta{eta_list[etaIND]}.dat','w' )
        for p in range( NPolCompute ):
            dipFile.write( f'{eff_dipole_1[ etaIND, p ]}\n')
        dipFile.close()
        #eff_dipole_1[etaIND,:] = np.loadtxt(f'data_dipole/dipole_Polaritons_eta{round(eta_list[etaIND],8)}_wc{round(wc,4)}.dat')

    return eff_dipole_1

def plot_TD( TD_matter, TD_pol ):

    NMatPlot = [ 6 ]
    NPolPlot = [ 6,7 ] #np.arange(NPolCompute)

    for j in NMatPlot: # Only 0-N transition density, Z-Dimension
        T = np.sum( TD_matter[0,j,:,:,:],axis=(1,2) )
        #plt.plot( np.arange(Nxyz[0]) * 0.529 * dLxyz[0], np.abs(T),'o', linewidth=6, alpha=0.5, label=f'$M_{j}^x$' )
        plt.plot( np.arange(Nxyz[0]) * 0.529 * dLxyz[0], np.abs(T),'-', linewidth=6, alpha=0.5, label=f'$M_{j}^x$' )
    for p in NPolPlot: # 0: P0 --> P0, Z-Dimensions
        for etaIND in range(Neta):
            eta = round(eta_list[etaIND],8)
            Px = np.sum( TD_pol[etaIND,p,:,:,:],axis=(1,2) )
            plt.plot( np.arange(Nxyz[0]) * 0.529 * dLxyz[0], np.abs(Px),'--', linewidth=2, label=f'$P0{p}^x (\eta = {eta})$' )
            
    plt.legend(loc='upper right')
    plt.xlabel('Length (A)',fontsize=15) 
    plt.ylabel('Transition Density',fontsize=15)
    plt.tight_layout()
    plt.savefig("data_TD/Transition_Density_Slices_POL_Matter_X.jpg")
    plt.clf()


    for j in NMatPlot: #Nad+1): # Only 0-N transition density
        T = np.sum( TD_matter[0,j,:,:,:],axis=(0,2) )
        #plt.plot( np.arange(Nxyz[1]) * 0.529 * dLxyz[1], np.abs(T),'o', linewidth=6, alpha=0.5, label=f'$M_{j}^y$' )
        plt.plot( np.arange(Nxyz[1]) * 0.529 * dLxyz[1], np.abs(T),'-', linewidth=6, alpha=0.5, label=f'$M_{j}^y$' )
    for p in NPolPlot: # 0: P0 --> P0
        for etaIND in range(Neta):
            eta = round(eta_list[etaIND],8)
            Py = np.sum( TD_pol[etaIND,p,:,:,:],axis=(0,2) )
            plt.plot( np.arange(Nxyz[1]) * 0.529 * dLxyz[1], np.abs(Py),'--', linewidth=2, label=f'$P0{p}^y (\eta = {eta})$' )
            
    plt.legend(loc='upper right')
    plt.xlabel('Length (A)',fontsize=15)
    plt.ylabel('Transition Density',fontsize=15)
    plt.tight_layout()
    plt.savefig("data_TD/Transition_Density_Slices_POL_Matter_Y.jpg")
    plt.clf()

    for j in NMatPlot: #Nad+1): # Only 0-N transition density
        T = np.sum( TD_matter[0,j,:,:,:],axis=(0,1) )
        #plt.plot( np.arange(Nxyz[2]) * 0.529 * dLxyz[2], np.abs(T),'o', linewidth=6, alpha=0.5, label=f'$M_{j}^z$' )
        plt.plot( np.arange(Nxyz[2]) * 0.529 * dLxyz[2], np.abs(T),'-', linewidth=6, alpha=0.5, label=f'$M_{j}^z$' )
    for p in NPolPlot: # 0: P0 --> P0
        for etaIND in range(Neta):
            eta = round(eta_list[etaIND],8)
            Pz = np.sum( TD_pol[etaIND,p,:,:,:],axis=(0,1) )
            plt.plot( np.arange(Nxyz[2]) * 0.529 * dLxyz[2], np.abs(Pz),'--', linewidth=2, label=f'$P0{p}^z (\eta = {eta})$' )
            
    plt.legend(loc='upper right')
    plt.xlabel('Length (A)',fontsize=15)
    plt.ylabel('Transition Density',fontsize=15)
    plt.tight_layout()
    plt.savefig("data_TD/Transition_Density_Slices_POL_Matter_Z.jpg")
    plt.clf()

def compute_Transition_density_1r( Upol, TD_matter ):
    TD_pol = get_TD_fast( Upol, TD_matter, NPolCompute ).reshape(( len(Upol), NPolCompute, Nxyz[0],Nxyz[1],Nxyz[2] ))

    if ( write_TD_Files == True ):
        # Print TD Files in Gaussian cube format
        for p in range(NPolCompute): # p = 1 --> T_{01}(R)
            for etaIND in range(len(eta_list)):
                eta = round(eta_list[etaIND],8)
                print (f"Writing transition density file. p: 0-->{p}, eta = {eta}")
                f = open(f'data_TD/trans_eta{eta}_Nad{Nad}_Nfock{NFock}_P0{p}.cube','w')
                f.write(f"P0{p+1} Transition Density\n")
                f.write(f"Totally {NGrid} grid points\n")
                f.write(f"{NAtoms} {Lxyz[0]} {Lxyz[1]} {Lxyz[2]}\n")
                f.write(f"{Nxyz[0]} {dLxyz[0]}  0.000000   0.000000\n")
                f.write(f"{Nxyz[1]} 0.000000   {dLxyz[1]} 0.000000\n")
                f.write(f"{Nxyz[2]} 0.000000   0.000000   {dLxyz[2]} \n")
                for at in range(len(coords)):
                    f.write( coords[at] )
                for x in range(Nxyz[0]):
                    #print(f'X = {x}')
                    for y in range(Nxyz[1]):
                        outArray = []
                        for z in range(Nxyz[2]):
                            outArray.append( TD_pol[etaIND,p,x,y,z] )
                            if ( len(outArray) % 6 == 0 or z == Nxyz[2]-1 ):
                                #outArray.append('\n')
                                f.write( " ".join(map( str, np.round(outArray,8) )) + "\n" )
                                outArray = []
                f.close()


    plot_TD( TD_matter, TD_pol )

def compute_Transition_density_1r1q( Upol, TD_matter ):
    dict_xyz = { 'x':0, 'y':1, 'z':2 }
    HO_WFNs = [ getHO_q(n) for n in range(NFock) ]
    TD_projected = getTraceTD( np.abs(TD_matter) ) # state, state, plotDIM
    TD_pol = get_TD_fast_1r1q( Upol, TD_projected, NPolCompute, HO_WFNs )

    RGrid = np.linspace( 0, Lxyz[ dict_xyz[plotDIM] ], Nxyz[ dict_xyz[plotDIM] ] )
    
    for etaIND in range( len(eta_list) ):
        print ( f' Saving contour plot for eta # {etaIND} ' )
        #for p in range( NPolCompute ):
        for p in [ 0, 1, 6, 7, 21, 22 ]: 
            plt.contourf( RGrid, QGrid, TD_pol[ etaIND, p, :, : ].T )
            plt.xlabel( r' r ($\AA$)', fontsize=15 )
            plt.ylabel( ' q$_c$ (a.u.)', fontsize=15 )
            plt.title( r' < P0 | $\hat\rho (r,q_c)$ | P'+f'{p} >', fontsize=15 )
            plt.colorbar()
            plt.tight_layout()
            plt.savefig( f'data_TD/TD_1r1q_eta{eta_list[etaIND]}_Nad{Nad}_Nfock{NFock}_P0{p}.jpg' )
            plt.clf()


def getDipole_pow_2(Upol,MU):       

    print("Starting dipole square matrix elements.")
    eff_dipole_2 = compute_dipole_2_fast( NFock, Upol, MU, NFockInclude, NMatterInclude, NPolCompute )
    np.savetxt(f'data_dipole/dipole_pow_2.dat', np.c_[ np.arange(len(eff_dipole_2)), eff_dipole_2 ] )
    return eff_dipole_2
        
def getOscStr(Epol,eff_dipole_1):        
        
    print("Starting oscillator strength calculations.")
    eff_osc_str = compute_eff_osc_str_fast( NPolar, Epol, eff_dipole_1 )
    np.savetxt(f'data_dipole/osc_str.dat', np.c_[ np.arange(len(eff_osc_str)), eff_osc_str ] )
    return eff_osc_str

def plotSpectra(Epol,eff_osc_str,Char):

    EMin = 0.5 # eV
    EMax = 3.0 # eV
    Npts = 1000

    sig = 0.1 # eV

    dE = (EMax - EMin) / Npts
    energy = np.linspace( EMin, EMax, Npts )

    Epol_transition = np.zeros(( Neta, NPolar ))
    for j in range( NPolar ):
        Epol_transition[:,j] = (Epol[:,j] - Epol[:,0])


    # Make SPECTRA with gaussian width sig
    Spec = np.zeros(( Neta, Npts ))
    for etaIND in range( len(eta_list) ):
        for k in range( len(energy) ):
            #Spec[etaIND,k] += np.sum( eff_osc_str[etaIND,:] * np.exp( -(energy[k] - Epol_transition[etaIND,:]) ** 2 / 2 / sig ** 2 ) )
            Spec[etaIND,k] += np.sum( eff_osc_str[etaIND,:] * (sig/4) * sig / ( (energy[k] - Epol_transition[etaIND,:])**2 + (0.5*sig)**2 ) )

    X, Y = np.meshgrid(energy, eta_list)
    #Spec[Spec < 0.0001] = -1 
    
    np.savetxt( f"data_spectra/polaritonic_SPECTRA_E{d}_wc{wc}_Nad{Nad+1}_Nfock{NFock}_Ndipoles{NPolCompute}_sig{sig}.dat", Spec )



    # GET FRANK STYLE

    """
    # Add white to beginning using space from first color
    rainbow = mpl.colormaps['rainbow']#.resampled(256)
    violet = rainbow( [0] )
    white = np.array([256/256, 256/256, 256/256, 1])
    NTotal = 10000
    NShift = 2000
    newcolors = rainbow(np.linspace(0, 1, NTotal))
    newcolors[NShift:,:] = rainbow(np.linspace(0, 1, NTotal - NShift))
    for j in range( NShift ):
        frac = j / NShift
        newcolors[j, :] = np.array( (1-frac) * white + frac * violet )
    cmap = ListedColormap(newcolors)
    """

    # Add white to beginning
    rainbow = mpl.colormaps['rainbow']#.resampled(256)
    violet = rainbow( [0] )
    white = np.array([256/256, 256/256, 256/256, 1])
    NTotal = 10000
    NShift = 100
    NStart = 50
    newcolors = rainbow(np.linspace(0, 1, NTotal))
    newcolors[NShift:,:] = rainbow(np.linspace(0, 1, NTotal - NShift))
    for j in range( NShift ):
        if( j > NStart ):
            frac = (j-NStart) / (NShift) * (NShift/NStart)
        else:
            frac = 0
        newcolors[j, :] = np.array( (1-frac) * white + frac * violet )
    cmap = ListedColormap(newcolors)




    #cmap = mpl.cm.magma
    #cmap = mpl.cm.magma_r
    #cmap = mpl.cm.Purples
    #cmap = mpl.cm.terrain_r
    #cmap = mpl.cm.rainbow
    #cmap = mpl.cm.PuBu
    #cmap = mpl.cm.afmhot_r
    #plt.contourf( X, Y, Spec, cmap=cmap, norm=mpl.colors.Normalize(vmin=0, vmax=35)  )
    #plt.contourf( X, Y, Spec, cmap=cmap )
    plt.imshow( Spec, origin='lower', interpolation='gaussian', extent=(np.amin(energy), np.amax(energy), np.amin(eta_list), np.amax(eta_list)), cmap=cmap, aspect='auto', norm=mpl.colors.Normalize(vmin=0, vmax=30) )
    #plt.imshow( Spec, origin='lower', extent=(np.amin(energy), np.amax(energy), np.amin(eta_list), np.amax(eta_list)), cmap=cmap, aspect='auto', norm=mpl.colors.Normalize(vmin=0, vmax=35) )

    print(np.shape(Spec))


    #plt.colorbar()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=30), cmap=cmap), pad=0.01)

    plt.xlim(0.5,3.0)
    #plt.ylim(eta_list[0],eta_list[-1])
    plt.xlabel("Energy (eV)", fontsize=10)
    plt.ylabel("eta",fontsize=10)
    plt.title(f"Absorption (NExc:{Nad} Nfock: {NFock})",fontsize=10)
    #plt.colorbar()
    plt.tight_layout()
    plt.savefig(f"data_spectra/polaritonic_SPECTRA_E{d}_wc{wc}_Nad{Nad+1}_Nfock{NFock}_Ndipoles{NPolCompute}.jpg", dpi=600)
    plt.clf()




    ############### Make DOS with gaussian width sig #################
    Npts = 1000
    sig = 0.01 # eV
    #sig_eta = 0.0001 # eV
    energy = np.linspace( EMin, EMax, Npts )
    #eta_grid = np.linspace( eta_list[0], eta_list[-1], Npts )
    DOS = np.zeros(( Neta, Npts ))
    #DOS = np.zeros(( Npts, Npts ))
    for k in range( len(energy) ):
        #print(f"{k+1} of = {len(energy)}" )   
        #for ek in range( len(eta_grid) ):         
        for etaIND, eta in enumerate( eta_list ):
            #DOS[etaIND,k] += np.sum( np.exp( -(energy[k] - Epol_transition[etaIND,:] ) ** 2 / 2 / sig ** 2 ) )
            #DOS[ek,k] += np.sum( np.exp( -(energy[k] - Epol_transition[etaIND,:] ) ** 2 / 2 / sig ** 2 ) * np.exp( -(eta_grid[ek] - eta_list[etaIND])**2 / 2 / sig_eta**2 ) )
            DOS[etaIND,k] += np.sum( np.sum( 1.0000 * (sig/4) * sig / ( (energy[k] - Epol_transition[etaIND,:])**2 + (0.5*sig)**2 ) ) )

    X, Y = np.meshgrid(energy, eta_list)
    #X, Y = np.meshgrid(energy, eta_grid)
    #DOS[DOS < 0.0001] = -1 
    



    #cmap = mpl.cm.magma
    #cmap = mpl.cm.Purples
    cmap = mpl.cm.terrain_r
    cmap = mpl.cm.rainbow
    #plt.contourf( X, Y, DOS, cmap=cmap )
    #plt.contourf( X, Y, DOS, cmap=cmap, norm=mpl.colors.Normalize(vmin=0, vmax=2)  )
    #plt.contourf( X, Y, np.abs(DOS), cmap=cmap, norm=mpl.colors.LogNorm()  )
    #plt.imshow( DOS , origin='lower')
    #plt.pcolormesh( X, Y, DOS, cmap=cmap, norm=mpl.colors.Normalize(vmin=0, vmax=4)  )
    plt.imshow( DOS, origin='lower', interpolation='gaussian', extent=(np.amin(energy), np.amax(energy), np.amin(eta_list), np.amax(eta_list)), cmap=cmap, aspect='auto', norm=mpl.colors.Normalize(vmin=0, vmax=1) )
    #plt.imshow( DOS, origin='lower', extent=(np.amin(energy), np.amax(energy), np.amin(eta_list), np.amax(eta_list)), cmap=cmap, aspect='auto', norm=mpl.colors.Normalize(vmin=0, vmax=1) )

    


    #plt.colorbar()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=1), cmap=cmap))
    #plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.LogNorm(vmin=0.1), cmap=cmap))

    plt.xlim(-0.05,EMax)
    plt.ylim(eta_list[0],eta_list[-1])
    plt.xlabel("Energy (eV)", fontsize=10)
    plt.ylabel("eta",fontsize=10)
    plt.title(f"DOS (NExc:{Nad} Nfock: {NFock})",fontsize=10)
    #plt.colorbar()
    plt.tight_layout()
    plt.savefig(f"data_spectra/polaritonic_DOS_E{d}_wc{wc}_Nad{Nad+1}_Nfock{NFock}_Ndipoles{NPolCompute}.jpg", dpi=800)
    plt.clf()

    # Produce character-colored plot of polariton energies
    # Let's just save the data nicely and plot in origin
    np.savetxt(f"data_spectra/eta_list.dat", eta_list[:] )
    np.savetxt(f"data_spectra/EPol_Transition.dat", Epol_transition[:,:] )
    np.savetxt(f"data_spectra/Photon_Number.dat", Char[:,:] )


    NINTERP = 5000
    NSTATE_INTERP = 25

    eta_interp = np.linspace( eta_list[0], eta_list[-1], NINTERP )
    Epol_interp = np.zeros(( NSTATE_INTERP, NINTERP ))
    Char_interp = np.zeros(( NSTATE_INTERP, NINTERP ))

    for state in range( NSTATE_INTERP ):
        f_x = interp1d( eta_list, Epol_transition[:,state] )
        Epol_interp[ state, : ] = f_x( eta_interp )

        f_x = interp1d( eta_list, Char[:,state] )
        Char_interp[ state, : ] = f_x( eta_interp )

    np.savetxt(f"data_spectra/eta_list_interp.dat", eta_interp[:] )
    np.savetxt(f"data_spectra/EPol_Transition_interp.dat", Epol_interp[:,:].T )
    np.savetxt(f"data_spectra/Photon_Number_interp.dat", Char_interp[:,:].T )




def main():
    getGlobals()
    Epol, Char, Upol, MU = get_HadMU()
    #TDM_matter = get_TD_Data()
    #TD_pol_1r    = compute_Transition_density_1r( Upol, TDM_matter ) 
    #TD_pol_1r1q  = compute_Transition_density_1r1q( Upol, TDM_matter ) 



    #getExansion(Upol)

    start_D1 = time()
    eff_dipole_1 = getDipole_pow_1(Upol,MU)
    end_D1 = time()
    #print (f'\t\tDipole_1: {round(end_D1-start_D1,2)} s')
    eff_osc_str = getOscStr(Epol,eff_dipole_1)
    plotSpectra(Epol,eff_osc_str,Char)



    #dipole_pow_2 = getDipole_pow_2(Upol,MU) # DONT NEED THIS
    
    
    
    #make_movie()


if ( __name__ == '__main__'):
    main()
