%chk=geometry.chk
%mem=20GB
%nprocshared=12

#p wB97XD/6-311++G** nosymm
#P TD=(singlets,nstates=500) IOp(6/8=3) IOp(9/40=6)

Title Card Required

0 1
C          0.04850        0.03916        0.15341
H         -0.88424       -0.51794        0.15089
C          1.23714       -0.63094        0.08168
C          2.49547        0.06104        0.08154
H          3.39864       -0.57505        0.02030
H          1.22913       -1.71091        0.02336
O          2.62769        1.28119        0.14326
N         -0.07781        1.37100        0.22960
H         -0.98063        1.80507        0.27940
H          0.77230        1.92511        0.23435







