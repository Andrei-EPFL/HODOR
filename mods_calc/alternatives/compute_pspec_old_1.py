import os
import configparser
import numpy as np

import MAS_library as MASL
import Pk_library as PKL

from nbodykit.source.catalog import ArrayCatalog
from nbodykit.lab import FFTPower

class ComputePKClass():
	""" The class that computes the pspec given the config file and a sample of galaxies """
	def __init__(self, config_file):
		config = configparser.ConfigParser()
		config.read(config_file)

		self.Lbox = config['params'].getfloat('Lbox')
		self.verbose = config['params'].getboolean('verbose')
		self.Nmesh = config['pspec'].getint('Nmesh')
		self.kmin = config['pspec'].getfloat('kmin')
		self.fit_kmin = config['pspec'].getfloat('fit_kmin')
		self.fit_kmax = config['pspec'].getfloat('fit_kmax')
		self.dk = config['pspec'].getfloat('dk')
		self.interlaced = config['pspec'].getboolean('interlaced')
		self.Nmu = config['pspec'].getint('n_mu_bins')

	def pspec_pylians(self, dict_of_gsamples, index):
		grid = self.Nmesh
		BoxSize = self.Lbox
		verbose = self.verbose
		MAS = "CIC"
		threads = 8
		axis = 2

		pk0_list = []
		pk2_list = []
		pk4_list = []
		
		for j, key in enumerate(dict_of_gsamples.keys()):
			sample1 = dict_of_gsamples[key]


			# define 3D density field
			delta = np.zeros((grid,grid,grid), dtype=np.float32)

			# construct 3D density field
			MASL.MA(sample1, delta, BoxSize, MAS, verbose=verbose)

			# at this point, delta contains the effective number of particles in each voxel
			# now compute overdensity and density constrast
			delta /= np.mean(delta, dtype=np.float64)
			delta -= 1.0

			# compute power spectrum
			Pk = PKL.Pk(delta, BoxSize, axis, MAS, threads, verbose)

			# 3D P(k)
			k       = Pk.k3D
			pk0_list.append(Pk.Pk[:,0])
			pk2_list.append(Pk.Pk[:,1])
			pk4_list.append(Pk.Pk[:,2])


		mean_pk0 = np.mean(np.array(pk0_list), 0)
		mean_pk2 = np.mean(np.array(pk2_list), 0)
		mean_pk4 = np.mean(np.array(pk4_list), 0)

		range_ = (self.fit_kmin <= k) & (k <= self.fit_kmax)
		return k[range_], mean_pk0[range_], mean_pk2[range_], mean_pk4[range_]


	def pspec_1d(self, dict_of_gsamples, index):
		pk_list = []
		
		for j, key in enumerate(dict_of_gsamples.keys()):
			sample1 = dict_of_gsamples[key]
			cat = ArrayCatalog(sample1)
			mesh = cat.to_mesh(Nmesh=self.Nmesh, interlaced=self.interlaced, BoxSize=self.Lbox, compensated=True)
			r = FFTPower(mesh, second=None, mode='1d', kmin=self.kmin, Nmu=self.Nmu, dk=self.dk, Nmesh=self.Nmesh)
			
			Pk = r.power
			pk = Pk['power'].real - Pk.attrs['shotnoise']
			
			pk_list.append(pk)

		mean_pk = np.mean(np.array(pk0_list), 0)

		return Pk['k'], mean_pk


	def pspec_2d(self, dict_of_gsamples, index):
		pk0_list = []
		pk2_list = []
		pk4_list = []
		
		for j, key in enumerate(dict_of_gsamples.keys()):
			sample1 = dict_of_gsamples[key]
			
			sample_xyz= np.zeros(len(sample1), dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4')])
			sample_xyz['x'] = sample1[:, 0]
			sample_xyz['y'] = sample1[:, 1]
			sample_xyz['z'] = sample1[:, 2]

			cat = ArrayCatalog(sample_xyz)
			cat['Position'] = cat['x'][:, None] * [1, 0, 0] + cat['y'][:, None] * [0, 1, 0] + cat['z'][:, None] * [0, 0, 1]

			mesh = cat.to_mesh(Nmesh=self.Nmesh, position="Position", interlaced=self.interlaced, BoxSize=self.Lbox, compensated=True)
			r = FFTPower(mesh, second=None, mode='2d', kmin=self.kmin, Nmu=self.Nmu, los=[0, 0, 1], poles=[0, 2, 4], dk=self.dk, Nmesh=self.Nmesh)

			poles = r.poles

			pk0 = poles['power_0'].real - poles.attrs['shotnoise']
			pk2 = poles['power_2'].real
			pk4 = poles['power_4'].real
			pk0_list.append(pk0)
			pk2_list.append(pk2)
			pk4_list.append(pk4)

		mean_pk0 = np.mean(np.array(pk0_list), 0)
		mean_pk2 = np.mean(np.array(pk2_list), 0)
		mean_pk4 = np.mean(np.array(pk4_list), 0)

		range_ = (self.fit_kmin <= poles['k']) & (poles['k'] <= self.fit_kmax)
		return poles['k'][range_], mean_pk0[range_], mean_pk2[range_], mean_pk4[range_]

