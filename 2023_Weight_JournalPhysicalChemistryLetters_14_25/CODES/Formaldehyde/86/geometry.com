%chk=geometry.chk
%mem=2GB
%nprocshared=1

#p B3LYP/6-311+G**
#p TD=(singlets,nstates=30)IOp(6/8=3) IOp(9/40=4)

Scan Step 86, Bond Length 0.9299999999999999

0 1
C	0.00000	0.00000	-0.52736	
H	0.00000	0.93726	-1.11871	
H	-0.00000	-0.93726	-1.11871	
O	-0.00000	-0.00000	0.4026399999999999	









