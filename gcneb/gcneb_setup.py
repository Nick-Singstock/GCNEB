#! /usr/bin/env python

'''
Script to create images for NEB jobs.

Author: Nick Singstock
'''

import argparse
from ase.io import read, write
from ase.neb import NEB
import os

# setup argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--initial',type=str, default='00/CONTCAR',
    help='The initial POSCAR or CONTCAR')
parser.add_argument('-f', '--final',type=str, default='NIMAGES/CONTCAR',
    help='The final POSCAR or CONTCAR')
parser.add_argument('-n', '--nimages',type=int, default=5,
    help='Number of images to make')
parser.add_argument('--idpp', help='Whether to use idpp interpolation (default True)',
                    default=True)
args = parser.parse_args()

# function for fixing POSCAR files
def insert_el(filename):
    with open(filename, 'r') as f:
        file = f.read()
    contents = file.split('\n')
    ele_line = contents[0]
    if contents[5].split() != ele_line.split():
        contents.insert(5, ele_line)
    with open(filename, 'w') as f:
        f.write('\n'.join(contents))

# get initial and final structures
st_init = read(args.initial, format='vasp')
del st_init.constraints

final = args.final
if final == 'NIMAGES/CONTCAR':
    final = str(args.nimages).zfill(2) + '/CONTCAR'
st_final = read(final, format='vasp') 

# add images for NEB as copies of initial state
images = [st_init.copy() for i in range(args.nimages + 1)]
images.append(st_final)

# interpolate new images
neb = NEB(images)
neb.interpolate(mic = True)
if args.idpp:
    neb.idpp_interpolate(traj = None, steps = 20, mic = True)

# save new images to folders
for i, image in enumerate(images):
    # skip initial and final folders 
    if i == 0 or i >= len(images)-1:
        continue
    folder = str(i).zfill(2)
    if not os.path.exists(folder):
        os.mkdir(folder)
    image.write(os.path.join(folder,'POSCAR'),format="vasp", direct=True)
    insert_el(os.path.join(folder,'POSCAR'))
   
print('\nGCNEB directory setup. Calculation can be run with sub_JDFTx.py and "nimages'+
      str(args.nimages) + '" in inputs file \n')
