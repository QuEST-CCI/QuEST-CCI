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
                        "H" ,\
                        "H" ,\
                        "H" ,\
                        "N" ,\
                        "O" ,\
                        "O" ,\
                        "C" ,\
                        "H" ,\
                        "C" ,\
                        "H" ,\
                        "Br"])
    COORDS = np.array([
                 [ 0.51932475,    1.23303451,   -0.03194925],
                 [ 1.94454413,    1.26916358,   -0.03672882],
                 [ 2.62037793,    0.09283428,   -0.02499003],
                 [-0.19603352,    0.03013062,    0.00102732],
                 [-0.02069420,    2.17423764,   -0.04336646],
                 [ 2.48281698,    2.20891057,   -0.03611879],
                 [-1.27770137,    0.03990295,    0.01166953],
                 [ 4.09213475,    0.09594076,    0.03662979],
                 [ 4.63930696,   -1.02169275,    0.14459220],
                 [ 4.66489883,    1.19839699,   -0.02327545],
                 [ 0.49428518,   -1.16712649,    0.02099746],
                 [-0.03251071,   -2.11492669,    0.05447935],
                 [ 1.96291176,   -1.21653219,   -0.02111314],
                 [ 2.44359113,   -1.96306433,    0.61513886],
                 [ 2.17304025,   -1.94912156,   -1.90618750]])
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
%mem=1GB\n\
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


