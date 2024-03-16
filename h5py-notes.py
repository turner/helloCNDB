h5py Notes

cndbf = h5py.File(name + '.cndb', 'w')

H = cndbf.create_group('Header')
H.attrs.update(info)

C = cndbf.create_group(info_dict['name'])
pos = C.create_group('spatial_position')

pos.create_dataset(str(frame[-1]+1), data=spatial_pos)


C.create_dataset('genomic_position',data=np.array(genseq))

C.create_dataset('time',data=np.array(frame))

H.attrs.update({'chromosome' : chr_name})

cndbf.close()
