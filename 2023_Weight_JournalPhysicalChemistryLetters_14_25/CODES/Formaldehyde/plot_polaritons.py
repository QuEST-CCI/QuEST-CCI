import numpy as np
from matplotlib import pyplot as plt
import sys


NSteps = 300 # 0, 1, 2, 3, ..., Nsteps;  # 40 = 5.0 Ang., 60 = 6.0 Ang.
NStatesPOL = 20#len( np.loadtxt( f"Epol_Ez_0.dat" ) )
NStatesAD = 3#len( np.loadtxt( f"Had_Ez_0.dat" ) )
Epol = np.zeros(( NSteps, NStatesPOL ))
photon_number = np.zeros(( NSteps, NStatesPOL ))
Ead = np.zeros(( NSteps, NStatesAD ))

A0 = round( float(sys.argv[1]) ,3)
wc  = round( float(sys.argv[2]) ,2)

data_dir = "pol_data"

#print ("Polaritonic States:", NStatesPOL)
#print ("Adiabatic States:", NStatesAD)

print ("\t\tPlotting.")

R = np.arange( 0.5, NSteps * 0.005 + 0.5, 0.005 )


for d in ['z']:

    for step in range( NSteps):
        #print (f"Step: {step}")
        Epol[step] = np.loadtxt( f"{data_dir}/Epol_E{d}_A0{np.round(A0,4)}_wc{np.round(wc,4)}_{step}.dat" )[:NStatesPOL] # In eV
        #photon_number[step] = np.loadtxt( f"{data_dir}/photon_number_E{d}_A0{np.round(A0,4)}_wc{np.round(wc,4)}_{step}.dat" )[:NStatesPOL] # In eV
        Ead[step] = np.loadtxt( f"{data_dir}/Had_E{d}_{step}.dat" )[:NStatesAD]

    EZero_Ad = np.min( Ead[:,0]) # This is ground state minimum
    EZero_Pol = EZero_Ad # Same zero for ad. and pol.
                
    #print (f"EZero_Ad = { EZero_Ad }")
    #print (f"EZero_Pol = { EZero_Pol }")

    
    #plt.plot( R, Epol[:,0] - EZero_Pol, ":", linewidth=6, dash_capstyle='round', label=f"Polaritonic E{d}" )
    #for state in range( 1, 10 ):
        #print (f"State: {state}")
    #    plt.plot( R, Epol[:,state] - EZero_Pol, "--", linewidth=6 )     # PLOT POLARITONIC

#plt.plot( R, Ead[:,0] - EZero_Ad, linewidth=6, alpha=0.5, label="Adiabatic" )
#for state in range( 1, NStatesAD ):
#    plt.plot( R, Ead[:,state] - EZero_Ad ) # PLOT ADIABATIC

#plt.xlabel("CO Bond Distance (A)", fontsize=15)
#plt.ylabel("Energy (eV)", fontsize=15)
####plt.xlim(1,1.8)
####plt.ylim(7.5,10) # meV

#plt.xlim(1.1,1.3)
#plt.ylim(7.5,8.5) # meV

#plt.legend()
#plt.savefig(f"polaritons_E{d}_A0{np.round(A0,4)}_wc{np.round(wc,4)}.jpg")
#plt.clf()


np.savetxt(f"image_data/Had.dat", Ead[:,:10] - EZero_Ad )
np.savetxt(f"image_data/polaritons__E{d}_A0{np.round(A0,4)}_wc{np.round(wc,4)}.dat", Epol[:,:10] - EZero_Pol )
#np.savetxt(f"image_data/photon_number_eta_{np.round(A0,4)}_wc_{np.round(wc,4)}.dat", photon_number[:,:10] )

# Save minimum differences between states between 1.3 < R(C-O) < 1.4 between states P2/P3 and P3/P4 (all of A1 symmetry).
indices = np.arange(160,181)
R = R[indices]
EDiff = np.zeros(( 2 ))
RMin = np.zeros(( 2 ))
#for state in range( 10 ):
EDiff[0] = np.min( Epol[indices,2] - Epol[indices,1] )
EDiff[1] = np.min( Epol[indices,3] - Epol[indices,2] )

RMin[0] = R[ np.argmin(Epol[indices,2] - Epol[indices,1]) ]
RMin[1] = R[ np.argmin(Epol[indices,3] - Epol[indices,2]) ]

np.savetxt(f"image_data/EPol_Diff__E{d}_A0{np.round(A0,4)}_wc{np.round(wc,4)}.dat", np.c_[RMin[:], 1000*EDiff[:]], fmt='%.5f' )





"""
### Plot Adiabatic States
R = np.arange( 0.5, NSteps * 0.005 + 0.5, 0.005 )
#plt.plot( R, Ead[:,0] - EZero_Ad, linewidth=1.5 )
for state in range( 1, NStatesAD ):
    plt.plot( R, (Ead[:,state] - EZero_Ad) / 1000, label=f"S{state}" ) # PLOT ADIABATIC
plt.xlabel("CO Bond Distance (A)", fontsize=15)
plt.ylabel("Energy (eV)", fontsize=15)
plt.legend()
plt.savefig("AdiabaticOnly.jpg")
plt.clf()
"""