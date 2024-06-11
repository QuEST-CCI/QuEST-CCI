import numpy as np
import subprocess as sp
from matplotlib import pyplot as plt
import os, sys


RPA = True # TDA is always printed... if RPA=True, look for RPA results
DATA_DIR = "PLOTS_DATA"

def plot_dipole(DIP,label):

    xyz_dict = { 0:'x', 1:'y', 2:'z' }
    for dim in range(3):
        plt.imshow( np.abs( DIP[:,:,dim] ),origin='lower', cmap='afmhot_r')
        plt.colorbar(pad=0.01)
        plt.xlabel("Electronic State Index",fontsize=15)
        plt.ylabel("Electronic State Index",fontsize=15)
        plt.title(f"Electronic Dipole Matrix, |$\mu^{xyz_dict[dim]}$"+"$_{\\alpha\\beta}$| (a.u.)",fontsize=15)
        plt.savefig(f"{DATA_DIR}/DIP_{label}_{xyz_dict[dim]}.jpg", dpi=600)
        plt.clf()

        np.savetxt( f"{DATA_DIR}/DIP_{label}_{xyz_dict[dim]}_0J.dat", DIP[0,:,dim]  )
        plt.stem( np.arange(len(DIP)), np.abs( DIP[0,:,dim] ), bottom=-10)
        plt.xlabel("Electronic State Index",fontsize=15)
        plt.ylabel("Dipole Squared, $(\hat{\mu}^2)_{0 \\beta}$",fontsize=15)
        plt.xlim(-0.5,len(DIP)+0.5)
        plt.ylim(-0.1)
        plt.savefig(f"{DATA_DIR}/DIP_{label}_{xyz_dict[dim]}_0J.jpg", dpi=600)
        plt.clf()

    plt.imshow( np.sqrt( DIP[:,:,0]**2 + DIP[:,:,1]**2 + DIP[:,:,2]**2 ),origin='lower', cmap='afmhot_r')
    plt.colorbar(pad=0.01)
    plt.xlabel("Electronic State Index",fontsize=15)
    plt.ylabel("Electronic State Index",fontsize=15)
    plt.title("Electronic Dipole Matrix, |$\mu^{Tot.}_{\\alpha\\beta}$| (a.u.)",fontsize=15)
    plt.savefig(f"{DATA_DIR}/DIP_{label}_Total.jpg", dpi=600)
    plt.clf()

    np.savetxt( f"{DATA_DIR}/DIP_{label}_Total_0J.dat", np.sqrt( DIP[0,:,0]**2 + DIP[0,:,1]**2 + DIP[0,:,2]**2  ) )
    plt.stem( np.arange(len(DIP)), np.sqrt( DIP[0,:,0]**2 + DIP[0,:,1]**2 + DIP[0,:,2]**2  ), bottom=-10)
    plt.xlabel("Electronic State Index",fontsize=15)
    plt.ylabel("Dipole, $\mu^{Tot.}_{0 \\beta}$",fontsize=15)
    plt.xlim(-0.5,len(DIP)+0.5)
    plt.ylim(-0.1)
    plt.savefig(f"{DATA_DIR}/DIP_{label}_Total_0J.jpg", dpi=600)
    plt.clf()

def plot_dipole_square(DIP,label):

    xyz_dict = { 0:'x', 1:'y', 2:'z' }
    DIP_SQUARE = np.einsum( "JKd,KLd->JLd" , DIP[:,:,:], DIP[:,:,:] )
    for dim in range(3):
        plt.imshow( np.abs( DIP_SQUARE[:,:,dim] ),origin='lower', cmap='afmhot_r')
        plt.colorbar(pad=0.01)
        plt.xlabel("Electronic State Index",fontsize=15)
        plt.ylabel("Electronic State Index",fontsize=15)
        plt.title("Dipole Matrix Squared, $(\hat{\mu}^2)_{\\alpha \\beta}$"+f"({xyz_dict[dim]}) (a.u.)",fontsize=15)
        plt.savefig(f"{DATA_DIR}/DIP_SQUARE_{label}_{xyz_dict[dim]}.jpg", dpi=600)
        plt.clf()

        np.savetxt( f"{DATA_DIR}/DIP_SQUARE_{label}_{xyz_dict[dim]}_0J.dat", DIP_SQUARE[0,:,dim]  )
        plt.stem( np.arange(len(DIP_SQUARE)), np.abs( DIP_SQUARE[0,:,dim] ), bottom=-10)
        plt.xlabel("Electronic State Index",fontsize=15)
        plt.ylabel("Dipole Squared, $(\hat{\mu}^2)_{0 \\beta}$",fontsize=15)
        plt.xlim(-0.5,len(DIP_SQUARE)+0.5)
        plt.ylim(-0.1)
        plt.savefig(f"{DATA_DIR}/DIP_SQUARE_{label}_{xyz_dict[dim]}_0J.jpg", dpi=600)
        plt.clf()



    plt.imshow( np.sqrt( DIP_SQUARE[:,:,0]**2 + DIP_SQUARE[:,:,1]**2 + DIP_SQUARE[:,:,2]**2 ),origin='lower', cmap='afmhot_r')
    plt.colorbar(pad=0.01)
    plt.xlabel("Electronic State Index",fontsize=15)
    plt.ylabel("Electronic State Index",fontsize=15)
    plt.title("Dipole Matrix Squared, $(\hat{\mu}^2)_{\\alpha \\beta}$"+f"(Total) (a.u.)",fontsize=15)
    plt.savefig(f"{DATA_DIR}/DIP_SQUARE_{label}_Total.jpg", dpi=600)
    plt.clf()

    np.savetxt( f"{DATA_DIR}/DIP_SQUARE_{label}_Total_0J.dat", np.sqrt( DIP_SQUARE[0,:,0]**2 + DIP_SQUARE[0,:,1]**2 + DIP_SQUARE[0,:,2]**2)  )
    plt.stem( np.arange(len(DIP_SQUARE)), np.sqrt( DIP_SQUARE[0,:,0]**2 + DIP_SQUARE[0,:,1]**2 + DIP_SQUARE[0,:,2]**2), bottom=-10)
    plt.xlabel("Electronic State Index",fontsize=15)
    plt.ylabel("Dipole Squared, $(\hat{\mu}^2)^{Tot.}_{0 \\beta}$",fontsize=15)
    plt.xlim(-0.5,len(DIP_SQUARE)+0.5)
    plt.ylim(-0.1)
    plt.savefig(f"{DATA_DIR}/DIP_SQUARE_{label}_Total_0J.jpg", dpi=600)
    plt.clf()

def plot_energies(E_TDA,E_RPA):
    dRPA = (E_RPA - E_RPA[0]) * 27.2114 # a.u. to eV
    dTDA = (E_TDA - E_TDA[0]) * 27.2114 # a.u. to eV
    plt.plot( np.arange(len(E_TDA)), dRPA, c="black", label="RPA" )
    plt.plot( np.arange(len(E_TDA)), dTDA, c="red", label="TDA" )
    plt.xlim(0,len(E_TDA))
    plt.xlabel("Electronic State Index",fontsize=15)
    plt.ylabel("Excitation Energy Difference, RPA - TDA (eV)",fontsize=15)
    plt.title("Electronic Dipole Matrix, |$\mu^{Tot.}_{\\alpha\\beta}$| (a.u.)",fontsize=15)
    plt.savefig(f"{DATA_DIR}/E_TDA_RPA_COMP.jpg", dpi=600)
    plt.clf()

    plt.plot( np.arange(len(E_TDA)), dRPA - dTDA )
    plt.xlabel("Electronic State Index",fontsize=15)
    plt.ylabel("Excitation Energy Difference, RPA - TDA (eV)",fontsize=15)
    plt.xlim(0,len(E_TDA))
    plt.title("Electronic Dipole Matrix, |$\mu^{Tot.}_{\\alpha\\beta}$| (a.u.)",fontsize=15)
    plt.savefig(f"{DATA_DIR}/E_TDA_RPA_DIFF.jpg", dpi=600)
    plt.clf()


def get_nuclear_dipole_old():
    # TODO -- Shift to center of charge before computing dipole for charges systems
    
    def get_atoms():

        with open("QCHEM.out","r") as file01:
            lines       = file01.readlines()
            atomFLAG    = False
            atom_coords = []
            atom_labels = []
            units       = ""
            for count, line in enumerate(lines):
                t = line.split()
                if ( len(t) > 1 and t[0] == "Standard" and t[1] == "Nuclear" ):
                    units = t[-1].split("(")[1].split(")")[0]
                    atomFLAG = True
                    counter = 3
                    while (atomFLAG == True):
                        S = lines[count+counter].split()
                        if ( len(S) == 5 ):
                            atom_labels.append( S[1] )
                            atom_coords.append( S[2:] )
                        if ( "----------------------------------------------------------------" in S ):
                            break
                        counter += 1

        return atom_labels, np.array(atom_coords).astype(float), units.upper()

    def get_nuclear_charges():
        dict = {
        "H":0,\
        "He":0,\
        "Li":0,\
        "Be":0,\
        "B":0,\
        "C":0,\
        "N":0,\
        "O":0,\
        "F":0,\
        "Ne":0,\
        "Na":0,\
        "Mg":0,\
        "Al":0,\
        "Si":0,\
        "P":0,\
        "S":0,\
        "Cl":0,\
        "Ar":0,\
        "K":0,\
        "Ca":0,\
        "Sc":0,\
        "Ti":0,\
        "V":0,\
        "Cr":0,\
        "Mn":0,\
        "Fe":0,\
        "Co":0,\
        "Ni":0,\
        "Cu":0,\
        "Zn":0,\
        "Ga":0,\
        "Ge":0,\
        "As":0,\
        "Se":0,\
        "Br":0,\
        "Kr":0,\
        "Rb":0,\
        "Sr":0,\
        "Y":0,\
        "Zr":0,\
        "Nb":0,\
        "Mo":0,\
        "Tc":0,\
        "Ru":0,\
        "Rh":0,\
        "Pd":0,\
        "Ag":0,\
        "Cd":0,\
        "In":0,\
        "Sn":0,\
        "Sb":0,\
        "Te":0,\
        "I":0,\
        "Xe":0,\
        "Cs":0,\
        "Ba":0,\
        "La":0,\
        "Ce":0,\
        "Pr":0,\
        "Nd":0,\
        "Pm":0,\
        "Sm":0,\
        "Eu":0,\
        "Gd":0,\
        "Tb":0,\
        "Dy":0,\
        "Ho":0,\
        "Er":0,\
        "Tm":0,\
        "Yb":0,\
        "Lu":0,\
        "Hf":0,\
        "Ta":0,\
        "W":0,\
        "Re":0,\
        "Os":0,\
        "Ir":0,\
        "Pt":0,\
        "Au":0,\
        "Hg":0,\
        "Tl":0,\
        "Pb":0,\
        "Bi":0,\
        "Po":0,\
        "At":0,\
        "Rn":0,}
        for count,at in enumerate(dict):
            dict[at] = float(count+1)

        return dict

    ANG_to_DEBYE = 4.80320
    DEBYE_to_AU  = 1/2.54158
    
    atom_labels, atom_coords, units = get_atoms()
    nuclear_charge_dict = get_nuclear_charges()

    # Get the center of nuclear charge and total charge
    COC    = np.zeros((3)) # Center of nuclear charge (hopefully close to overall center)
    Q_T    = 0.0 # Total nuclear charge
    for at, atom in enumerate( atom_labels ):
        COC[0] += nuclear_charge_dict[atom] * atom_coords[at,0]
        COC[1] += nuclear_charge_dict[atom] * atom_coords[at,1]
        COC[2] += nuclear_charge_dict[atom] * atom_coords[at,2]
        Q_T    += nuclear_charge_dict[atom]

    # Get the arbitrary origin nuclear dipole moment
    MU_NUC = np.zeros((3))
    for at, atom in enumerate( atom_labels ):
        MU_NUC[0] += nuclear_charge_dict[atom] * atom_coords[at,0]
        MU_NUC[1] += nuclear_charge_dict[atom] * atom_coords[at,1]
        MU_NUC[2] += nuclear_charge_dict[atom] * atom_coords[at,2]

    # Test the nuclear dipole shift to nuclear center of charge
    #print("\nTESTING NUCLEAR CHARGE SHIFT TO COC:")
    #print( "Nuclear Dipole:\n\t", MU_NUC*ANG_to_DEBYE * DEBYE_to_AU )
    #print( "Nuclear Dipole (Shifted to COC):\n\t", (MU_NUC - COC/Q_T)*ANG_to_DEBYE * DEBYE_to_AU )
    #print( "\n")

    if ( units == "Angstroms".upper() ):
        MU_NUC *= ANG_to_DEBYE * DEBYE_to_AU
    elif ( units == "Bohrs".upper() ):
        print("Only works for Angstrom*|e| for now.")
        exit()
    else:
        print("Only works for Angstrom*|e| for now.")
        exit()
    
    #print("Nuclear Dipole:", MU_NUC)

    return MU_NUC

def get_nuclear_dipole( MU_00_EL ):

    DEBYE_to_AU  = 1/2.54158
    MU_NUC       = np.zeros((3))
    
    # Get dX,dY,dZ of total (el + nuc) ground state dipole from QCHEM output
    MU_00_TOT    = np.array(sp.check_output(''' grep "Dipole Moment (Debye)" QCHEM.out -A 1 | tail -n 1 | awk '{print $2, $4, $6}' ''', shell=True).decode().split()).astype(float)
    MU_00_TOT   *= DEBYE_to_AU
    MU_NUC       = MU_00_TOT - MU_00_EL

    print("OLD (a.u.) :", MU_NUC)
    print("NEW (a.u.) :", -get_nuclear_dipole_old())
    print("DIFF (a.u.):", MU_NUC + get_nuclear_dipole_old())

    return MU_NUC

sp.call(f"mkdir -p {DATA_DIR}", shell=True)

##### READ QCHEM OUTPUT #####

# Look for ground state total SCF energy
GS_ENERGY = float( sp.check_output("grep 'SCF   energy in the final basis set' QCHEM.out | awk '{print $9}'", shell=True).decode().split()[0] )

# Look for excitation energies from TDA and RPA
sp.call("grep 'excitation energy' QCHEM.out > "+f"{DATA_DIR}/EXC_ENERGIES.dat", shell=True)
NROOTS = int( sp.check_output(f"wc -l {DATA_DIR}/EXC_ENERGIES.dat "+ "| awk '{print $1}'", shell=True).decode() )

if ( RPA == True ): 
    NROOTS //= 2
    print( f"\tThere are {NROOTS} TDDFT/TDA roots." )
    print( f"\tThere are {NROOTS} TDDFT/RPA roots." )
    sp.call(f"head -n {NROOTS} {DATA_DIR}/EXC_ENERGIES.dat" + " | awk '{print $8}' > "+f"{DATA_DIR}/EXC_TDA.dat", shell=True)
    sp.call(f"tail -n {NROOTS} {DATA_DIR}/EXC_ENERGIES.dat" + " | awk '{print $8}' > "+f"{DATA_DIR}/EXC_RPA.dat", shell=True)
else:
    print( f"\tThere are {NROOTS} TDDFT/TDA roots." )
    sp.call(f"head -n {NROOTS} {DATA_DIR}/EXC_ENERGIES.dat" + " | awk '{print $8}' > "+f"{DATA_DIR}/EXC_TDA.dat", shell=True)
sp.call(f"rm {DATA_DIR}/EXC_ENERGIES.dat",shell=True)

##### Create arrays to store all data #####
EAD     = np.zeros(( NROOTS+1 )) # Add ground state
DIP     = np.zeros(( NROOTS+1, NROOTS+1, 3 )) # Add ground state
#MU_NUC  = get_nuclear_dipole_old()
###########################################
EAD[0] = GS_ENERGY # Already in a.u.
if ( RPA == True ):
    EAD[1:] = GS_ENERGY + np.loadtxt(f"{DATA_DIR}/EXC_TDA.dat") / 27.2114 # eV to a.u.
    np.savetxt(f"{DATA_DIR}/ADIABATIC_ENERGIES_TDA.dat",EAD)
    E_TDA = EAD * 1.0
    EAD[1:] = GS_ENERGY + np.loadtxt(f"{DATA_DIR}/EXC_RPA.dat") / 27.2114 # eV to a.u.
    np.savetxt(f"{DATA_DIR}/ADIABATIC_ENERGIES_RPA.dat",EAD)
    E_RPA = EAD * 1.0
    plot_energies(E_TDA,E_RPA)
else:
    EAD[1:] = GS_ENERGY + np.loadtxt(f"{DATA_DIR}/EXC_TDA.dat") / 27.2114 # eV to a.u.
    np.savetxt(f"{DATA_DIR}/ADIABATIC_ENERGIES_TDA.dat",EAD)


##### Look for dipoles from TDA #####
DIP[0,0,:] = np.array( sp.check_output("grep -A 4 'Electron Dipole Moments of Ground State' QCHEM.out | tail -n 1", shell=True).decode().split()[1:] ).astype(float)
DIP[np.arange(1,NROOTS+1), np.arange(1,NROOTS+1),:]  = np.array( sp.check_output(f"grep -m 1 -A {3+NROOTS} 'Electron Dipole Moments of Singlet Excited State' QCHEM.out | tail -n {NROOTS}", shell=True).decode().split() ).astype(float).reshape(NROOTS,4)[:,1:]
DIP[0, np.arange(1,NROOTS+1),:] = np.array( sp.check_output(f"grep -m 1 -A {3+NROOTS} 'Transition Moments Between Ground and Singlet Excited States' QCHEM.out | tail -n {NROOTS}", shell=True).decode().split() ).astype(float).reshape(NROOTS,6)[:,2:5]
DIP[np.arange(1,NROOTS+1),0,:] = np.array( sp.check_output(f"grep -m 1 -A {3+NROOTS} 'Transition Moments Between Ground and Singlet Excited States' QCHEM.out | tail -n {NROOTS}", shell=True).decode().split() ).astype(float).reshape(NROOTS,6)[:,2:5]

# Excited-to-excited needs more work to be done...need a loop
NTERMS = NROOTS * (NROOTS-1) // 2 # Upper triangle without diagonal terms
tmp = np.array( sp.check_output(f"grep -m 1 -A {3+NTERMS} 'Transition Moments Between Singlet Excited States' QCHEM.out | tail -n {NTERMS}", shell=True).decode().split() ).astype(float).reshape(NTERMS,6)[:,:5]
for J,K,dx,dy,dz in tmp:
    DIP[int(J),int(K),:] = np.array([dx,dy,dz])
    DIP[int(K),int(J),:] = np.array([dx,dy,dz])

# Calculate nuclear dipole from QCHEM GS total (el + nuc) dipole
MU_NUC  = get_nuclear_dipole( DIP[0,0,:] )

# Add nuclear dipole (which is already of opposite sign)
for state in range(NROOTS+1):
    DIP[state,state,:] += MU_NUC

np.save(f"{DATA_DIR}/DIPOLE_TDA.dat.npy", DIP)
plot_dipole(DIP,"TDA")
plot_dipole_square(DIP,"TDA")

if ( RPA == True ):

    DIP = np.zeros(( NROOTS+1, NROOTS+1, 3 )) # Add ground state

    ##### Look for dipoles from RPA #####
    DIP[0,0,:] = np.array( sp.check_output("grep -A 4 'Electron Dipole Moments of Ground State' QCHEM.out | tail -n 1", shell=True).decode().split()[1:] ).astype(float)
    DIP[np.arange(1,NROOTS+1), np.arange(1,NROOTS+1),:] = np.array( sp.check_output(f"grep -m 2 -A {3+NROOTS} 'Electron Dipole Moments of Singlet Excited State' QCHEM.out | tail -n {NROOTS}", shell=True).decode().split() ).astype(float).reshape(NROOTS,4)[:,1:]
    DIP[0, np.arange(1,NROOTS+1),:] = np.array( sp.check_output(f"grep -m 2 -A {3+NROOTS} 'Transition Moments Between Ground and Singlet Excited States' QCHEM.out | tail -n {NROOTS}", shell=True).decode().split() ).astype(float).reshape(NROOTS,6)[:,2:5]
    DIP[np.arange(1,NROOTS+1),0,:] = np.array( sp.check_output(f"grep -m 2 -A {3+NROOTS} 'Transition Moments Between Ground and Singlet Excited States' QCHEM.out | tail -n {NROOTS}", shell=True).decode().split() ).astype(float).reshape(NROOTS,6)[:,2:5]

    # Excited-to-excited needs more work to be done...need a loop
    NTERMS = NROOTS * (NROOTS-1) // 2 # Upper triangle without diagonal terms
    tmp = np.array( sp.check_output(f"grep -m 2 -A {3+NTERMS} 'Transition Moments Between Singlet Excited States' QCHEM.out | tail -n {NTERMS}", shell=True).decode().split() ).astype(float).reshape(NTERMS,6)[:,:5]
    for J,K,dx,dy,dz in tmp:
        DIP[int(J),int(K),:] = np.array([dx,dy,dz])
        DIP[int(K),int(J),:] = np.array([dx,dy,dz])

    # Add nuclear dipole (with negative sign) from perm. dipoles:
    # Nuclear dipole is same for RPA and TDA (maybe obviously...)
    for state in range(NROOTS+1):
        DIP[state,state,:] += MU_NUC
    
    np.save(f"{DATA_DIR}/DIPOLE_RPA.dat.npy", DIP)
    np.savetxt(f"{DATA_DIR}/DIPOLE_RPA_x.dat", DIP[:,:,0])
    plot_dipole(DIP,"RPA")
    plot_dipole_square(DIP,"RPA")





