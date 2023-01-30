#!/usr/bin/env python3
import numpy as np
import bigfile
import h5py
from nbodykit.transform import HaloConcentration, HaloRadius
from nbodykit.cosmology import Planck15
import sys

if len(sys.argv) != 3:
  print('Usage: {} INPUT OUTPUT'.format(sys.argv[0]), file=sys.stderr)
  sys.exit(1)


ifile = sys.argv[1]
ofile = sys.argv[2]

with bigfile.File(ifile) as f:
  # data = bigfile.Dataset(f['LL-0.200/'], ['ID','Length','Position','Velocity'])
  data = bigfile.Dataset(f['LL-0.200/'], ['Length','Position','Velocity'])
  header = f['Header']
  header2 = f['LL-0.200']


redshift = 1.0 / header.attrs['ScalingFactor'][0] - 1
particle_mass = header2.attrs['M0'][0] * 1e10
print("particle mass = ", particle_mass)
cosmo = Planck15.match(Omega0_m=header.attrs['OmegaM'])

mass = data['Length'][:] * particle_mass
print("min number of particles = ", np.min(data['Length'][:]))
print("max number of particles = ", np.max(data['Length'][:]))
Conc = HaloConcentration(mass, cosmo, redshift)
Rvir = HaloRadius(mass, cosmo, redshift)
Rs = Rvir / Conc
upid = -np.ones(data.size)

# import collections
# import json
# counter=collections.Counter(data['Length'][:])
# x=[]
# y=[]
# for key in counter.keys():
#   y.append(int(counter[key]))
#   x.append(int(key))

# x = np.array(x)
# y = np.array(y)
# np.savetxt(ofile, np.array([x, y]).T)


with h5py.File(ofile, 'w') as ff:
  ff.create_group('halo')
  # ff.create_dataset('halo/ID', data=data['ID'][:])
  ff.create_dataset('halo/ID', data=np.arange(len(mass)))
  ff.create_dataset('halo/Mvir', data=mass, dtype=np.float32)
  ff.create_dataset('halo/Rvir', data=Rvir, dtype=np.float32)
  ff.create_dataset('halo/Rs', data=Rs, dtype=np.float32)
  ff.create_dataset('halo/X', data=data['Position'][:][:,0])
  ff.create_dataset('halo/Y', data=data['Position'][:][:,1])
  ff.create_dataset('halo/Z', data=data['Position'][:][:,2])
  ff.create_dataset('halo/VX', data=data['Velocity'][:][:,0])
  ff.create_dataset('halo/VY', data=data['Velocity'][:][:,1])
  ff.create_dataset('halo/VZ', data=data['Velocity'][:][:,2])
  ff.create_dataset('halo/PID', data=upid, dtype=np.int64)
