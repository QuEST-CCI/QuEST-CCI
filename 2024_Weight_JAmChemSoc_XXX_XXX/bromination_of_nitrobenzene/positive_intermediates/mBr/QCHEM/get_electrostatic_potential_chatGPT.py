import numpy as np
import argparse

# Define the parser for command line arguments
parser = argparse.ArgumentParser(description='Compute the electrostatic potential from a Gaussian cube file.')
parser.add_argument('filename', type=str, help='The name of the Gaussian cube file.')
parser.add_argument('--output', type=str, default='potential.cube', help='The name of the output cube file.')
args = parser.parse_args()

# Read in the Gaussian cube file
with open(args.filename, 'r') as f:
    lines = f.readlines()

# Extract the relevant information from the Gaussian cube file
natoms = int(lines[2].split()[0])
origin = np.array([x for x in lines[2].split()[1:4]]).astype(float)
npoints = np.array([ lines[3].split()[0], lines[4].split()[0], lines[5].split()[0] ]).astype(int)
spacing = np.array([lines[3].split()[1], lines[4].split()[2], lines[5].split()[3]]).astype(float)
density = np.zeros(npoints)
for i in range(npoints[0]):
    for j in range(npoints[1]):
        for k in range(npoints[2]):
            index = 6 + natoms + i*npoints[1]*npoints[2] + j*npoints[2] + k
            try:
                t = lines[index].split()
                for item in t:
                    density[i,j,k] = float(item)
            except IndexError:
                continue

# Calculate the charge density from the electron density
charge_density = -1.0/(4.0*np.pi) * np.gradient(np.gradient(np.gradient(density, spacing[2], axis=2), spacing[1], axis=1), spacing[0], axis=0)

# Compute the electrostatic potential from the charge density
potential = np.zeros(npoints)
for i in range(npoints[0]):
    for j in range(npoints[1]):
        for k in range(npoints[2]):
            potential[i,j,k] = np.sum(charge_density[:i+1,:j+1,:k+1]) * np.product(spacing)

# Output the electrostatic potential to a new Gaussian cube file
with open(args.output, 'w') as f:
    f.write(lines[0])
    f.write(lines[1])
    f.write(lines[2])
    f.write('{} {} {} {}\n'.format(natoms, origin[0], origin[1], origin[2]))
    f.write('{} {} {} {}\n'.format(npoints[0], spacing[0], 0.0, 0.0))
    f.write('{} {} {} {}\n'.format(npoints[1], 0.0, spacing[1], 0.0))
    f.write('{} {} {} {}\n'.format(npoints[2], 0.0, 0.0, spacing[2]))
    for i in range(npoints[0]):
        for j in range(npoints[1]):
            outArray = []
            for k in range(npoints[2]):
                #f.write('{}\n'.format(potential[i,j,k]))
                outArray.append( potential[i,j,k] )
                if ( len(outArray) % 6 == 0 or k == npoints[2]-1 ):
                    f.write( " ".join(map( str, np.round(outArray,8)*1e6 )) + "\n" )
                    outArray = []