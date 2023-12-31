import numpy as np
import imageio
import subprocess as sp
from pygifsicle import optimize as gifOPT # This needs to be installed somewhere
from PIL import Image, ImageDraw, ImageFont

"""
Install pygifcicle:
pip3 install pygifsicle


Install gifsicle: ( echo "$(pwd)" = /scratch/bweight/software/ )
curl -sL http://www.lcdf.org/gifsicle/gifsicle-1.91.tar.gz | tar -zx
cd gifsicle-1.91
./configure --disable-gifview
make install exec_prefix=$(pwd) prefix=$(pwd) datarootdir=$(pwd)
"""

FontFile="/scratch/bweight/anaconda3/lib/python3.8/site-packages/anaconda_navigator/static/fonts/UbuntuMono-R.ttf"

def addText(fileloc,text):
    img = Image.open(fileloc)
    d1 = ImageDraw.Draw(img)
    d1.text( (400, 150), text, (0,0,0), font=ImageFont.truetype(FontFile,25) )
    img.save(f"{filename.split('.')[0]}_withTEXT.jpg")

#eta_list = np.round( np.arange( 0.005, 0.035, 0.01 ) ,3)
wc_list  = np.round( np.arange( 6.50, 10.01, 0.01  ) ,2)

movieNAME = 'movie.gif'
with imageio.get_writer(movieNAME, mode='I', fps=10) as writer: # Get a writer object
    #for eta in eta_list:
    eta = 0.04
    for wc in wc_list:
        filename = f"polaritons_eta_{eta}_wc_{wc}.jpg"
        print ("Filename:", filename)
        #continue
        try:
            addText(filename,f"A0 = {eta} a.u.\nwc = {wc} eV")
            image = imageio.imread( f"{filename.split('.')[0]}_withTEXT.jpg" ) # Read JPEG file
            sp.call( f" rm {filename.split('.')[0]}_withTEXT.jpg ",shell=True)
        except:
            print (f"Problem with {filename}. Skipping.")
            continue
        writer.append_data(image) # Write JPEG file (to memory at first; then printed at end)

gifOPT(movieNAME) # This will compress the GIF movie by at least a factor of two/three. With this: ~750 frames --> 80 MB

