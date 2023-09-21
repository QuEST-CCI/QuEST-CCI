%oldchk="../geometry.chk"
%chk=S_16.chk
%mem=10GB
%nprocshared=8

#p B3LYP/chkbasis
#p Geom=AllCheck Guess=(Read,Only) Density=(Check,Transition=16) Pop=(Minimal,NTO,SaveNTO)

TitleMe

0 1




