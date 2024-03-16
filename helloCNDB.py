#!/usr/bin/env python
#coding: utf8 

__description__ = \
"""
Converting *.sw to *.cndb format
"""

__author__ = " Antonio B. Oliveira Junior"
__date__   = "Sep/2023"

################################################################
# 
# Trajectories file *.sw to Compacted Nucleome Data Bank format .cndb
#
# usage:
#  ./ndb2cndb.py -f file.ndb -n name_CNDB_file
#
################################################################

import time
import argparse
import numpy as np
import h5py
from datetime import datetime

parser = argparse.ArgumentParser(description='Converting *.sw to *.cndb format')
parser.add_argument('-f', metavar='input-file-ndb-frames',help='sw file with frames',type=argparse.FileType('rt'))
parser.add_argument('-n', action='store', default='chromatin', dest='arg_name',  help='Name of output file')

try:
    arguments = parser.parse_args()
    print('################################################')
    print('Chosen file: {:}'.format(arguments.f.name))

except IOError as msg:
    parser.error(str(msg))                    

def to_float(value):
    try:
        return float(value)
    except ValueError:  # This will catch cases where the conversion fails, e.g., if value is 'nan'
        return float(0)

b_time = time.time()


file_sw = arguments.f
name    = arguments.arg_name
cndbf = h5py.File(name + '.cndb', 'w')


info_dict = {}
first_line = file_sw.readline().strip()
   
if first_line.startswith('#'):
    entries = first_line[2:].split()  #[2:] because have ##
    
    for entry in entries:
        key, value = entry.split('=')
        info_dict[key] = value

info = {
    'version' : '1.0.0',
    'info' : 'Encode',
    'title': 'The Nucleome Data Bank: Web-based Resources Simulate and Analyze the Three-Dimensional Genome',
    'expdta' : '',
    'author' : 'Antonio B Oliveira Junior',
    'cycle' : '',
    'date' : str(datetime.now()),
    'chains' : ''
    }

info.update(info_dict) #including information from the sw header

H = cndbf.create_group('Header')
H.attrs.update(info)

print('Converting file...')

for line in file_sw:
    info = line.split()

    if info[0] == 'chromosome':
        loop = 0
        types = []
        frame= [] 
        chr_name = ''
        inframe = False
        C = cndbf.create_group(info_dict['name'])
        #C.create_dataset('loops', data=[]) #empy parameter because its not prersent in sw file
        #C.create_dataset('types',data=[]) #empy parameter because its not prersent in sw file
        pos = C.create_group('spatial_position')

    elif info[0] == 'trace':
        if (inframe):
            print(frame[-1]+1)
            spatial_pos = np.column_stack((x, y, z))
            pos.create_dataset(str(frame[-1]+1), data=spatial_pos)
            inframe = False

        frame.append(int(info[1]))
        geni, genj = [], []
        x,y,z = [], [], []
        n_atoms = 0
        inframe = True

    elif len(info) == 6:
        geni.append(info[1])
        genj.append(info[2])
        x.append(to_float(info[3]))
        y.append(to_float(info[4]))
        z.append(to_float(info[5]))
        chr_name = info[0]
        print(chr_name)
        n_atoms += 1


spatial_pos = np.column_stack((x, y, z))
pos.create_dataset(str(frame[-1]+1), data=spatial_pos)

genseq = [[int(geni[n]), int(genj[n])] for n in range(len(geni))]
C.create_dataset('genomic_position',data=np.array(genseq))
        
C.create_dataset('time',data=np.array(frame))

H.attrs.update({'chromosome' : chr_name})

cndbf.close()

try:
    import hdf5_indexer
    print("indexing...")        
    hdf5_indexer.make_index(name + '.cndb')
    print("cndb indexed")
except ValueError:
    print("Not indexed, missing hdf5_indexer")

print('Finished!')

e_time = time.time()
elapsed = e_time - b_time
print('Ran in %.3f sec' % elapsed)
