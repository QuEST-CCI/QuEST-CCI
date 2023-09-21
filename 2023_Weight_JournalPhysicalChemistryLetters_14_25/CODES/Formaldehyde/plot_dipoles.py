import numpy as np
from matplotlib import pyplot as plt

NSteps = 300
NStates = 30 # From calculation
NInclude = 13 # To plot
include = np.zeros(( NSteps,NStates+1 )) # Filter by symmetry

DIP = np.zeros(( NSteps,NStates+1,NStates+1,3 ))
AD = np.zeros(( NSteps,NStates+1 ))

for step in range(NSteps):
    print (f"Step = {step}")

    include[step,0] = 1 # Always include ground state
    for count, line in enumerate( open(f"{step}/all_Exc_Sym.dat","r").readlines() ): 
        sym = line.split('-')[1].split('\n')[0]
        if ( sym == 'A1' ): include[step,count+1] = 1
        #include[step,count+1] = 1 # Turn all on

    #print (include[step,:NInclude])
    for line in np.loadtxt(f"{step}/dipole_matrix_E.dat"):
        i = int( line[0] )
        j = int( line[1] )
        if ( include[step,i] == 1 and include[step,j] == 1 ):
            DIP[step,i,j,:] = line[2:5] 
            DIP[step,j,i,:] = DIP[step,i,j,:]
            #print( "\tIncluded:", i,j )
        else:
            DIP[step,i,j,:] = line[2:5] * np.float('NaN')
            DIP[step,j,i,:] = DIP[step,i,j,:]

    for line in np.loadtxt(f"{step}/adiabtic_energies.dat"):
        i = int( line[0] )
        AD[step,i] = float( line[1] )

"""
# Remove 0.0 values from plot
DIP_ORIG = DIP * 1.0
for step in range(NSteps):
    for j in range(NStates):
        for k in range(NStates):
            for d in range(3):
                if ( DIP[step,j,k,d] == 0 ):
                    DIP[step,j,k,d] = np.nan
"""

R = np.linspace( 0.5, NSteps * 0.005 + 0.5, NSteps )
AD = AD - np.min( AD[:,0] )
for d in ['x','y','z']:
    pol = np.array([ d=='x',d=='y',d=='z' ]) * 1.0
    for statej in range( NInclude ):
        for statek in range( statej, NInclude ):
            #if ( np.max(np.abs(DIP[:,statej,statek,np.argmax(pol)])) > 0.01 ):

            if ( statej <= NInclude + 1 and statek <= NInclude + 1 ):

                if ( statej == 0 and statek != 0 ): # g --> e, First-Order Interaction
                    plt.plot( R[:], np.abs(DIP[:,statej,statek,np.argmax(pol)]), "--", linewidth=3 )
                elif( statej == statek ): # Sj --> Sj, Perm. Dipoles
                    plt.plot( R[:], np.abs(DIP[:,statej,statek,np.argmax(pol)]), "-", linewidth=2 )
                elif( statej != statek ): #  ej --> ek, 2nd-order Interaction
                    plt.plot( R[:], np.abs(DIP[:,statej,statek,np.argmax(pol)]), "-.", linewidth=2 )


    plt.xlabel("CO Bond Distance (A)", fontsize=15)
    plt.ylabel(f"|d{d}|", fontsize=15)
    plt.ylim(0,1.5)
    plt.xlim(1.28,1.35)  #(1,1.8)
    plt.legend()
    plt.savefig(f"DIP_vs_R_E{d}.jpg")
    plt.clf()

"""
AD_ROT = AD * 0.0 #  Copy shape of original matrix
for step in range( NSteps ):
    Edup, Udip = np.linalg.eigh( DIP_ORIG[ step, :,:,2 ] )
    print ( step, np.shape( AD[step,:] ), np.shape( Udip )  )
    AD_ROT[step,:] = AD[step,:] @ Udip.T

for state in range(2):
    #plt.plot( R, AD[:,state],"-",label=f"S{state}" )
    plt.plot( R, AD_ROT[:,state],"--",label=f"Dip{state}" )
plt.legend()
plt.ylim(-20,20)
plt.xlim(1,1.8)
plt.savefig("AD_vs_R.jpg")
"""