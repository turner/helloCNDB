#!/usr/bin/env python
#coding: utf8 

__description__ = \
"""
Converting *.sw to *.cndb format
"""

__author__ = " Antonio B. Oliveira Junior and Douglass Turner"
__date__   = "Mar/2024"

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


spacewalkFile = arguments.f
name    = arguments.arg_name
cndbf = h5py.File(name + '.cndb', 'w')


spacewalkMetaData = {}
first_line = spacewalkFile.readline().strip()
   
if first_line.startswith('#'):
    entries = first_line[2:].split()  #[2:] because have ##
    
    for entry in entries:
        key, value = entry.split('=')
        spacewalkMetaData[key] = value

metaData = {
    'version' : '1.0.0',
    'info' : 'Encode',
    'title': 'The Nucleome Data Bank: Web-based Resources Simulate and Analyze the Three-Dimensional Genome',
    'expdta' : '',
    'author' : 'Antonio B Oliveira Junior',
    'cycle' : '',
    'date' : str(datetime.now()),
    'chains' : ''
    }

metaData.update(spacewalkMetaData)

header = cndbf.create_group('Header')
header.attrs.update(metaData)

print('Converting file...')

for line in spacewalkFile:
    record = line.split()

    if record[0] == 'chromosome':
        loop = 0
        types = []
        frame= []
        inframe = False

        rootName = spacewalkMetaData['name']
        root = cndbf.create_group(rootName)

        genomicPosition = root.create_group('genomic_position')
        spatialPosition = root.create_group('spatial_position')

    elif record[0] == 'trace':
        if (inframe):
            print(frame[-1] + 1)

            genomicStack = np.column_stack((bpStart, bpEnd))
            genomicPosition.create_dataset(str(frame[-1] + 1), data=genomicStack)

            xyzStack = np.column_stack((x, y, z))
            spatialPosition.create_dataset(str(frame[-1] + 1), data=xyzStack)

            inframe = False

        frame.append(int(record[1]))
        bpStart, bpEnd = [], []
        x,y,z = [], [], []
        traceLength = 0
        inframe = True

    elif len(record) == 6:

        # genomic stack
        bpStart.append(int(record[1]))
        bpEnd.append(int(record[2]))

        # spatial stack
        x.append(to_float(record[3]))
        y.append(to_float(record[4]))
        z.append(to_float(record[5]))

        traceLength += 1

genomicStack = np.column_stack((bpStart, bpEnd))
genomicPosition.create_dataset(str(frame[-1] + 1), data=genomicStack)

xyzStack = np.column_stack((x, y, z))
spatialPosition.create_dataset(str(frame[-1] + 1), data=xyzStack)

# TODO: Add this to the Header attributes?
bpStartLength = len(bpStart)
rangeBPStartLength = range(bpStartLength)
genomicExtentList = [[int(bpStart[n]), int(bpEnd[n])] for n in range(len(bpStart))]

root.create_dataset('time', data=np.array(frame))

# TODO: Store chromosome in header?
# header.attrs.update({'chromosome' : chr})

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
