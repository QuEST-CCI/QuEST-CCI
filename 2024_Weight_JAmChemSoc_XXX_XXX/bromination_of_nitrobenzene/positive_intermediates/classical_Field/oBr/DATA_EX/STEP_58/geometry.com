%chk=geometry.chk
%mem=1GB
%nprocshared=12

#P wB97XD/6-311+G* NoSymm
#P Field=X+272

TitleMe

1 1
C 0.5193 1.2330 -0.0319
C 1.9445 1.2692 -0.0367
C 2.6204 0.0928 -0.0250
C -0.1960 0.0301 0.0010
H -0.0207 2.1742 -0.0434
H 2.4828 2.2089 -0.0361
H -1.2777 0.0399 0.0117
N 4.0921 0.0959 0.0366
O 4.6393 -1.0217 0.1446
O 4.6649 1.1984 -0.0233
C 0.4943 -1.1671 0.0210
H -0.0325 -2.1149 0.0545
C 1.9629 -1.2165 -0.0211
H 2.4436 -1.9631 0.6151
Br 2.1730 -1.9491 -1.9062






