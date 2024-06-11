import numpy as np
from glob import glob
from scipy.interpolate import interp2d

dA = 0.2
THETA_LIST = np.arange( 0, np.pi+dA, dA )
PHI_LIST   = np.arange( 0, np.pi+dA, dA )

A0 = 0.3
WC = 1.8






##### Compute oBr - mBr optimal theta/phi #####
E = np.zeros( (2, len(THETA_LIST), len(PHI_LIST)) ) # oBr, mBr
for THETAi,THETA in enumerate( THETA_LIST ):
    #print(f"Reading {THETAi+1} of {len(THETA_LIST)}")
    for PHIi,PHI in enumerate( PHI_LIST ):
        THETA = round(THETA,2)
        PHI   = round(PHI,2)
        E[0,THETAi, PHIi] = np.loadtxt( f"oBr/QCHEM/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{A0}_WC_{WC}_NF_5_NM_50.dat" )[0] / 27.2114 * 630
        E[1,THETAi, PHIi] = np.loadtxt( f"mBr/QCHEM/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{A0}_WC_{WC}_NF_5_NM_50.dat" )[0] / 27.2114 * 630

dE = E[0,:,:] - E[1,:,:]
# CRITICAL_ANGLE = np.unravel_index( np.argmin(dE), dE.shape )
# print("Critical Point (oBr/mBr):")
# print( "\tTheta = %1.1f (deg.)" % (THETA_LIST[CRITICAL_ANGLE[0]] * 180/np.pi) )
# print( "\tPhi   = %1.1f (deg.)" % (PHI_LIST[CRITICAL_ANGLE[1]] * 180/np.pi) )

dE_func    = interp2d( THETA_LIST, PHI_LIST, dE, kind='cubic' )
THETA_FINE = np.linspace( THETA_LIST[0], THETA_LIST[-1], 2000 )
PHI_FINE   = np.linspace( PHI_LIST[0], PHI_LIST[-1], 2000 )
dE_FINE    = dE_func( THETA_FINE, PHI_FINE ) 
CRITICAL_ANGLE = np.unravel_index( np.argmin(dE_FINE), dE_FINE.shape )

print("Critical Point (oBr/mBr):")
print( "\tTheta = %1.1f (deg.)" % (THETA_FINE[CRITICAL_ANGLE[0]] * 180/np.pi) )
print( "\tPhi   = %1.1f (deg.)" % (PHI_FINE[CRITICAL_ANGLE[1]] * 180/np.pi) )







##### Compute pBr - mBr optimal theta/phi #####
E = np.zeros( (2, len(THETA_LIST), len(PHI_LIST)) ) # oBr, mBr
for THETAi,THETA in enumerate( THETA_LIST ):
    #print(f"Reading {THETAi+1} of {len(THETA_LIST)}")
    for PHIi,PHI in enumerate( PHI_LIST ):
        THETA = round(THETA,2)
        PHI   = round(PHI,2)
        E[0,THETAi, PHIi] = np.loadtxt( f"pBr/QCHEM/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{A0}_WC_{WC}_NF_5_NM_50.dat" )[0] / 27.2114 * 630
        E[1,THETAi, PHIi] = np.loadtxt( f"mBr/QCHEM/PF_XYZ_SCAN/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{A0}_WC_{WC}_NF_5_NM_50.dat" )[0] / 27.2114 * 630

dE = E[0,:,:] - E[1,:,:]
# CRITICAL_ANGLE = np.unravel_index( np.argmin(dE), dE.shape )
# print("Critical Point (oBr/mBr):")
# print( "\tTheta = %1.1f (deg.)" % (THETA_LIST[CRITICAL_ANGLE[0]] * 180/np.pi) )
# print( "\tPhi   = %1.1f (deg.)" % (PHI_LIST[CRITICAL_ANGLE[1]] * 180/np.pi) )

dE_func    = interp2d( THETA_LIST, PHI_LIST, dE, kind='cubic' )
THETA_FINE = np.linspace( THETA_LIST[0], THETA_LIST[-1], 2000 )
PHI_FINE   = np.linspace( PHI_LIST[0], PHI_LIST[-1], 2000 )
dE_FINE    = dE_func( THETA_FINE, PHI_FINE ) 
CRITICAL_ANGLE = np.unravel_index( np.argmin(dE_FINE), dE_FINE.shape )

print("Critical Point (oBr/mBr):")
print( "\tTheta = %1.1f (deg.)" % (THETA_FINE[CRITICAL_ANGLE[0]] * 180/np.pi) )
print( "\tPhi   = %1.1f (deg.)" % (PHI_FINE[CRITICAL_ANGLE[1]] * 180/np.pi) )
