import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import os
import imageio
from pygifsicle import optimize as gifOPT # This needs to be installed somewhere
from PIL import Image, ImageDraw, ImageFont
import subprocess as sp

"""
Install pygifcicle:
pip3 install pygifsicle


Install gifsicle: ( echo "$(pwd)" = /scratch/bweight/software/ )
curl -sL http://www.lcdf.org/gifsicle/gifsicle-1.91.tar.gz | tar -zx
cd gifsicle-1.91
./configure --disable-gifview
make install exec_prefix=$(pwd) prefix=$(pwd) datarootdir=$(pwd)
"""

PIC_DIR = './IMAGES'
DATA_DIR = './pol_data'

eta_list = np.arange(0.0,0.04,0.001) # Coupling Strength
wc_list  = np.arange(1.0,2.2,0.02) # Photon Energy (eV)

NStates = 100

Epol = np.zeros(( len(eta_list), len(wc_list), NStates ))
Phot_Char = np.zeros(( len(eta_list), len(wc_list), NStates )) # Photonic character for each state

Neta = len(eta_list)
Nwc = len(wc_list)

for d in ['z']:
    for j in range(Neta):
        eta = np.round( eta_list[j],3 )
        print ( f"\tWorking on eta = {eta}" )
        for k in range(Nwc):
            wc = np.round( wc_list[k],2 )
            #print( f"\tWorking on wc = {wc}" )
            Epol[j,k] = np.loadtxt(f"{DATA_DIR}/Epol_E{d}_eta{eta}_wc{wc}.dat")[:NStates]
            Phot_Char[j,k,:] = np.loadtxt(f"{DATA_DIR}/Char_E{d}_eta{eta}_wc{wc}.dat")[:NStates]
            """
            with open( f"{DATA_DIR}/Upol_E{d}_eta{eta}_wc{wc}.dat","r") as f:
                for count,line in enumerate(f):
                    if ( count == 1 ):
                        Phot_Char[j,k,:] = (np.array(line.split()).astype(float) ** 2)[:NStates]
                        break
            """
if ( not os.path.exists(PIC_DIR) ):
    os.mkdir(PIC_DIR)

# Make Movie of Eta Scan
MOVE_NAME = 'ETA_SCAN.gif'
def makeJPG( filename, etaIND ):
    fig, ax = plt.subplots()

    for j in range( NStates ):
        im = ax.scatter( wc_list, Epol[ etaIND,:,j ],c=Phot_Char[etaIND,:,j],s=70,marker="o",edgecolor='none',cmap=plt.cm.cividis)
    eta = np.round( eta_list[etaIND],3 )
    plt.title(f"eta = {eta}")
    plt.xlim( wc_list[0],wc_list[-1] )
    plt.ylim(-0.05,3.0)
    plt.xlabel('Cavity Energy wc (eV)', fontsize=15)
    plt.ylabel('Polaritonic Energy (eV)', fontsize=15)
    fig.colorbar(im, ax=ax)
    im.set_clim(0.0, 1.0)
    fig.savefig(f'{filename}')
    plt.clf()

with imageio.get_writer(f"{PIC_DIR}/{MOVE_NAME}", mode='I', fps=2) as writer: # Get a writer object
    for etaIND in range( Neta ):
        eta = np.round( eta_list[etaIND],3 )
        filename = f'{PIC_DIR}/polariton_eta{eta}_wcSCAN.jpg'
        makeJPG( filename, etaIND )
        image = imageio.imread( f"{filename}" ) # Read JPEG file
        writer.append_data(image) # Write JPEG file (to memory at first; then printed at end)

sp.call(f"~/gifsicle-1.91/bin/gifsicle -b -O2 {PIC_DIR}/{MOVE_NAME}", shell=True)
#gifOPT(f"{PIC_DIR}/{MOVE_NAME}") # This will compress the GIF movie by at least a factor of two/three. With this: ~750 frames --> 80 MB

MOVE_NAME = 'WC_SCAN.gif'
def makeJPG( filename, wcIND ):
    fig, ax = plt.subplots()

    for j in range( NStates ):
        im = ax.scatter( eta_list * 1000, Epol[ :,wcIND,j ],c=Phot_Char[:,wcIND,j],s=70,marker="o",edgecolor='none',cmap=plt.cm.cividis)
    wc = np.round( wc_list[wcIND],3 )
    plt.title(f"wc = {wc} eV")
    plt.xlim( eta_list[0] * 1000 ,eta_list[-1] * 1000 )
    plt.ylim(-0.05,3.0)
    plt.xlabel('Coupling Strength eta (x10**3)', fontsize=15)
    plt.ylabel('Polaritonic Energy (eV)', fontsize=15)
    fig.colorbar(im, ax=ax)
    im.set_clim(0.0, 1.0)
    fig.savefig(f'{filename}')
    plt.clf()

with imageio.get_writer(f"{PIC_DIR}/{MOVE_NAME}", mode='I', fps=2) as writer: # Get a writer object
    for wcIND in range( Nwc ):
        wc = np.round( wc_list[wcIND],3 )
        filename = f'{PIC_DIR}/polariton_wc{wc}_etaSCAN.jpg'
        makeJPG( filename, wcIND )
        image = imageio.imread( f"{filename}" ) # Read JPEG file
        writer.append_data(image) # Write JPEG file (to memory at first; then printed at end)

sp.call(f"~/gifsicle-1.91/bin/gifsicle -b -O2 {PIC_DIR}/{MOVE_NAME}", shell=True)
#gifOPT(f"{PIC_DIR}/{MOVE_NAME}") # This will compress the GIF movie by at least a factor of two/three. With this: ~750 frames --> 80 MB




