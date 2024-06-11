import numpy as np
import subprocess as sp
import os

SUBMIT = True



def read_XYZ():

    LABELS = np.array([\
                        "C",\
                        "C" ,\
                        "C" ,\
                        "C" ,\
                        "C" ,\
                        "H" ,\
                        "H" ,\
                        "H" ,\
                        "H" ,\
                        "N" ,\
                        "O" ,\
                        "O" ,\
                        "C" ,\
                        "H" ,\
                        "Br"])
    COORDS = np.array([
                 [ 0.02949981,    1.33972592,    0.06817723],
                 [ 1.43483278,    1.28667967,    0.00635313],
                 [ 2.11179024,    0.05106117,   -0.00544138],
                 [ 1.44506636,   -1.13720058,    0.03116583],
                 [-0.68793171,    0.16822220,    0.10995314],
                 [-0.47126997,    2.29839666,    0.07811355],
                 [ 2.02732783,    2.19651728,   -0.03220624],
                 [ 1.98966526,   -2.07643217,    0.02318494],
                 [-1.77163480,    0.18040547,    0.15819632],
                 [ 3.58635895,    0.05097292,   -0.06745286],
                 [ 4.14711759,   -1.05966097,   -0.08807849],
                 [ 4.14497859,    1.16390951,   -0.09010823],
                 [-0.02361177,   -1.14582791,    0.08353483],
                 [-0.43674996,   -1.87247364,    0.78889576],
                 [-0.53591638,   -1.86972195,   -1.74078671]])
    return LABELS, COORDS

def make_COM( LABELS,COORDS, chi, dim ):

    if ( np.sign(chi) < 0 ):
        SIGN = "-"
    else:
        SIGN = "+"

    NATOMS, _ = COORDS.shape
    
    FILE01 = open("geometry.com","w")
    FILE01.write( 
f"\
%chk=geometry.chk\n\
%mem=9GB\n\
%nprocshared=12\n\n\
#P wB97XD/6-311+G* NoSymm\n\
#P Field={dim}{SIGN}{abs(int(chi))}\n\n\
TitleMe\n\n\
1 1\n" )

    for at in range( NATOMS ):
        FILE01.write( "%s %1.4f %1.4f %1.4f\n" % (LABELS[at], COORDS[at,0], COORDS[at,1], COORDS[at,2]) )
    FILE01.write("\n\n\n\n\n\n")


def make_directories( LABELS, COORDS, CHI  ):

    for dim in ["X","Y","Z"]:

        sp.call( f"mkdir -p DATA_E{dim}", shell=True )
        os.chdir(f"DATA_E{dim}/")

        NATOMS, _ = COORDS.shape
        for step,chi in enumerate(CHI):
            sp.call( f"mkdir -p STEP_{step}/", shell=True )
            os.chdir(f"STEP_{step}/")
            make_COM( LABELS, COORDS[:,:], chi, dim )


            sp.call( f"cp ../../submit.gaussian .", shell=True )
            if ( SUBMIT == True ):
                print ( f"Submitting step {step}" )
                sp.call( f"sbatch submit.gaussian", shell=True )
            os.chdir(f"../")
            
        os.chdir(f"../")






def main():
    LABELS, COORDS = read_XYZ()
    CHI = np.arange(-15,15,0.5) / 514 * 10000 # V/nm --> a.u. * 10,000
    make_directories(  LABELS, COORDS, CHI )

if ( __name__ == "__main__" ):
    main()


