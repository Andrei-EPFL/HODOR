import os
import sys
import configparser
# import multiprocessing as mp
from copy import copy
import numpy as np
#from functools import partial
import h5py
from halotools.sim_manager import UserSuppliedHaloCatalog
from halotools.empirical_models import HodModelFactory
from halotools.empirical_models import TrivialPhaseSpace, NFWPhaseSpace
from .hod_models import MWCens, MWSats, ContrerasCens, ContrerasSats

class ModelClass():
	""" The class that populates the halo catalog with galaxies
	given the config file """
	def __init__(self, config_file, halo_files):
		config = configparser.ConfigParser()
		config.read(config_file)

		self.model = config['params'].getint('model')
		if self.model == 0:
			mname = 'MW'
		elif self.model == 1:
			mname = 'GP18'
		else:
			raise ValueError('Unrecognised model')
		self.parameters_names = tuple(map(str, config.get('params', 'modelparamsnames').split(', ')))
		self.npar = config['params'].getint('npar')

		print('INFO: The used HOD model is {} having {} as parameters'.format(mname, self.parameters_names))
		if len(self.parameters_names) != self.npar:
			print("ERROR: Check the npar and modelparamsnames in the config file")
			sys.exit(1)

		self.verbose = config['params'].getboolean('verbose')
		self.redshift = config['params'].getfloat('redshift')
		self.Lbox = config['params'].getfloat('Lbox')
		self.Om = config['params'].getfloat('Om')
		self.initseed = config['params'].getint('initseed')
		self.Nseeds = config['params'].getint('Nseeds')

		print('\nINFO: The number of seeds to populate galaxies is: ', self.Nseeds)
		self.particle_mass = config['params'].getfloat('particle_mass')
		self.zspace = config['params'].getboolean('zspace')
		self.halo_files = halo_files
		print('\nINFO: The number of halo catalogs used to populate galaxies is: ', len(halo_files))
		
		self.halo_cat_list = self.compute_halocat()
		self.model_instance = self.compute_model_instance()
		#self.param0 = self.compute_param0()
		self.rsd_shift = self.compute_rsd_shift()
		self.Num_ptcl_requirement = config['params'].getfloat('Num_ptcl_requirement')

	def compute_halocat(self):
		""" Initialize a halo catalog """
		halo_cat_list = []
		
		for halo_file in self.halo_files:
			
			print(f'INFO: Compute the halo catalog: {os.path.basename(halo_file)}')
			f = h5py.File(halo_file, 'r')
			data = f['halo']
			print('INFO: min and max of X: ', np.min(data['X'][()]), np.max(data['X'][()]))
			print('INFO: min and max of Y: ', np.min(data['Y'][()]), np.max(data['Y'][()]))
			print('INFO: min and max of Z: ', np.min(data['Z'][()]), np.max(data['Z'][()]))
			X_data = data['Z'][()]
			Y_data = data['Y'][()]
			Z_data = data['X'][()]
			# Making all X,Y,Z values within 0-Lbox
			X_data[X_data < 0] = X_data[X_data < 0] + self.Lbox
			X_data[X_data > self.Lbox] = X_data[X_data > self.Lbox] - self.Lbox
			Y_data[Y_data < 0] = Y_data[Y_data < 0] + self.Lbox
			Y_data[Y_data > self.Lbox] = Y_data[Y_data > self.Lbox] - self.Lbox
			Z_data[Z_data < 0] = Z_data[Z_data < 0] + self.Lbox
			Z_data[Z_data > self.Lbox] = Z_data[Z_data > self.Lbox] - self.Lbox
			print('INFO: min and max of X: ', np.min(X_data), np.max(X_data))
			print('INFO: min and max of Y: ', np.min(Y_data), np.max(Y_data))
			print('INFO: min and max of Z: ', np.min(Z_data), np.max(Z_data))
			# Initialize the halo catalog
			halocat = UserSuppliedHaloCatalog(redshift=self.redshift, Lbox=self.Lbox, \
				particle_mass=self.particle_mass, halo_upid=data['PID'][()], \
				halo_x=X_data, halo_y=Y_data, halo_z=Z_data, \
				halo_vx=data['VZ'][()], halo_vy=data['VY'][()], halo_vz=data['VX'][()], \
				halo_id=data['ID'][()], halo_rvir=data['Rvir'][()], \
				halo_rs=data['Rs'][()], halo_mvir=data['Mvir'][()], \
				halo_nfw_conc=data['Rvir'][()]/data['Rs'][()], \
				halo_hostid=data['ID'][()])
			print('INFO: The length of the data is (i.e. the number of halos in the input file is): {:d}'.format(len(data['X'])))
			f.close()

			halo_cat_list.append(halocat)
			del data
		
		return halo_cat_list

	def compute_model_instance(self):
		""" Compute model instance """
		print('STATUS: Constructing HOD model ...')
		if self.model == 0:
			cens_occ_model = MWCens(redshift=self.redshift)
			sats_occ_model = MWSats(redshift=self.redshift)
		elif self.model == 1:
			cens_occ_model = ContrerasCens(redshift=self.redshift)
			sats_occ_model = ContrerasSats(redshift=self.redshift)

		cens_prof_model = TrivialPhaseSpace(redshift=self.redshift)
		sats_prof_model = NFWPhaseSpace(redshift=self.redshift)
		sats_occ_model._suppress_repeated_param_warning = True

		model_instance = HodModelFactory( \
			centrals_occupation=cens_occ_model, \
			centrals_profile=cens_prof_model, \
			satellites_occupation=sats_occ_model, \
			satellites_profile=sats_prof_model)
		return model_instance

	def compute_param0(self):
		""" Initialize default parameters """
		param0 = [self.model_instance.param_dict[self.parameters_names[i]] for i in range(self.npar)]
		print('INFO: The default parameters: ', param0)

		return param0

	def compute_rsd_shift(self):
		""" Compute the RSD redshift """
		if self.zspace:
			Ol = 1.0 - self.Om
			hubble = 100 * np.sqrt(self.Om * (self.redshift + 1.0)**3 + Ol)
			rsd_shift = (self.redshift + 1.0) / hubble
		else:
			rsd_shift = 0
		return rsd_shift

	def populate_mock(self, mfac, ref_num_dens):
		""" Populate the halo catalog with galaxies for different seeds """
		if self.verbose:
			print('\nSTATUS: Populate the halo catalog with galaxies...')
		for i in range(self.npar):
			self.model_instance.param_dict[self.parameters_names[i]] = mfac[i]

		# print('INFO: The values of the parameters {} are {}'.format(self.parameters_names, mfac))
		#if self.verbose:
		#	print('INFO: The new parameters: ', [self.model_instance.param_dict[self.parameters_names[i]] for i in range(self.npar)])
		seed0 = self.initseed
		Nseeds = self.Nseeds

		### Function which tells you how many central and satellite galaxies were alocated to the FastPM halo catalog.
		#self.check_galaxy_types(seed0)
		#exit()

		### Actually populating the halo mocks with galaxies for different seeds.
		return_dict = {}
		return_dict_gal_type = {}
		for j, halo_cat in enumerate(self.halo_cat_list):
			for i in range(Nseeds):
				seedt, sample1t, gal_type = self.intermediate_populate_mock(seed0 + i*100 + j, halo_cat)
				# seedt, sample1t, gal_type = self.intermediate_populate_mock(seed0 + i*1000 + indx, halo_cat)

				return_dict[os.path.basename(self.halo_files[j]) + str(seedt)] = sample1t
				return_dict_gal_type[seedt] = gal_type

				if (i == 0) and (j == 0):
					meas_num_dens = len(return_dict[os.path.basename(self.halo_files[j]) + str(seedt)]) / (self.Lbox ** 3)
					if np.abs(meas_num_dens - ref_num_dens) > 0.20 * ref_num_dens: ### Check whether one catalog is 8 Sigma_1 larger than the reference. Sigma_1 is the error of one catalog. 
						return return_dict, return_dict_gal_type

		return return_dict, return_dict_gal_type

	def intermediate_populate_mock(self, seed, halo_cat):
		""" This is an intermediate step in populating the halo catalog with galaxies,
		required to parallelize the procedure for multiple seeds """
		if self.verbose:
			print('SubSTATUS: Populate the halo catalog with galaxies using {} as seed...'.format(seed))
		
		self.model_instance.populate_mock(halo_cat, Num_ptcl_requirement=self.Num_ptcl_requirement, seed=seed)
		gtable = self.model_instance.mock.galaxy_table
		idx = (gtable['gal_type'] == 'centrals')
		gtable['gal_type'] = idx.astype(int)

		### Velocity Dispersion
		range_ = gtable['gal_type'] == 0
		vdisp = self.model_instance.param_dict["vdisp"]
		gtable['vz'][range_] = (gtable['vz'][range_] - gtable['halo_vz'][range_]) * vdisp + gtable['halo_vz'][range_]

		# RSD shift along OZ axis.
		gtable['z'] += gtable['vz'] * self.rsd_shift
		gtable['z'] %= self.Lbox

		sample1 = np.vstack((gtable['x'], gtable['y'], gtable['z'])).T
		if self.verbose:
			print('INFO: The numpy shape of the galaxy sample is {}'.format(sample1.shape))
		
		return seed, sample1, gtable['gal_type']

	def check_galaxy_types(self, seed):
		""" Function which tells you how many central and satellite galaxies were alocated to the FastPM halo catalog. """
		if self.verbose:
			print('SubSTATUS: Populate the halo catalog with galaxies using {} as seed...'.format(seed))
		
		for halo_cat in self.halo_cat_list:

			self.model_instance.populate_mock(halo_cat, Num_ptcl_requirement=self.Num_ptcl_requirement, seed=seed)
			gtable = self.model_instance.mock.galaxy_table

			idx = (gtable['gal_type'] == 'centrals')
			gtable['gal_type'] = idx.astype(int)

			unique, frequency = np.unique(gtable['gal_type'], return_counts=True)
			# print unique values array
			print("Unique Values:", unique)

			# print frequency array
			print("Frequency Values:", frequency)

	def return_galaxy_types(self, dict_gal_type):
		""" Function which tells you how many central and satellite galaxies were alocated to the FastPM halo catalog. """
		if self.verbose:
			print('\nSTATUS: Compute the 2pcf...')

		nsat = np.zeros(len(dict_gal_type.keys()))
		ncen = np.zeros(len(dict_gal_type.keys()))

		for j, key in enumerate(dict_gal_type.keys()):
			if self.verbose:
				print('INFO: >> Seed = {}...'.format(key))

			gal_type = dict_gal_type[key]

			unique, frequency = np.unique(gal_type, return_counts=True)
			# print unique values array
			#print("Unique Values:", unique)
			nsat[j] = frequency[0]

			if len(frequency) == 2:
				ncen[j] = frequency[1]
			# print frequency array
			#print("Frequency Values:", frequency)

		return np.mean(nsat), np.mean(ncen)

	# def parallel_populate_mock(self, mfac, ref_num_dens):
	# 	""" Populate the halo catalog with galaxies for different seeds in PARALLEL
	# 	Check populate_mock() function for eventual modifications!
	# 	"""
	# 	if self.verbose:
	# 		print('\nSTATUS: Populate in parallel the halo catalog with galaxies...')

	# 	for i in range(self.npar):
	# 		self.model_instance.param_dict[self.parameters_names[i]] = mfac[i]

	# 	# print('INFO: The values of the parameters {} are {}'.format(self.parameters_names, mfac))
	# 	#if self.verbose:
	# 	#	print('INFO: The new parameters: ', [self.model_instance.param_dict[self.parameters_names[i]] for i in range(self.npar)])

	# 	seed0 = self.initseed
	# 	Nseeds = self.Nseeds
	# 	if self.verbose:
	# 		print("INFO: There are {} available cores.".format(mp.cpu_count()))
	# 	if (mp.cpu_count() < Nseeds):
	# 		print("WARNING: The asked number of seeds is larger then the max number of cores so it is fixed to the max number of cores.")
	# 		Nseeds = mp.cpu_count()

	# 	manager = mp.Manager()
	# 	return_dict = manager.dict()
	# 	return_dict_gal_type = manager.dict()
	# 	jobs = []

	# 	### Actually populating the halo mocks with galaxies for different seeds.
	# 	for j, halo_cat in enumerate(self.halo_cat_list):
	# 		for i in range(Nseeds):
	# 			if (i == 0) and (j == 0):
	# 				self.parallel_intermediate_populate_mock(seed0 + i*100 + j, halo_cat, self.halo_files[j], return_dict, return_dict_gal_type)
	# 				meas_num_dens = len(return_dict[os.path.basename(self.halo_files[j]) + str(seed0 + i*100 + j)]) / (self.Lbox ** 3)
	# 				print("meas_dens: ", meas_num_dens)
	# 				if meas_num_dens > 1.2 * ref_num_dens:
	# 					return return_dict, return_dict_gal_type

	# 			else:	
	# 				p = mp.Process(target=self.parallel_intermediate_populate_mock, args=(seed0 + i*100 + j, halo_cat, self.halo_files[j], return_dict, return_dict_gal_type))
	# 				jobs.append(p)
	# 				p.start()

	# 	for proc in jobs:
	# 		proc.join()

	# 	return return_dict, return_dict_gal_type

	# def parallel_intermediate_populate_mock(self, seed, halo_cat, halo_name, return_dict, return_dict_gal_type):
	# 	""" This is an intermediate step in populating the halo catalog with galaxies,
	# 	required to parallelize the procedure for multiple seeds """
	# 	if self.verbose:
	# 		print('SubSTATUS: Populate the halo catalog with galaxies using {} as seed...'.format(seed))

	# 	self.model_instance.populate_mock(halo_cat, Num_ptcl_requirement=self.Num_ptcl_requirement, seed=seed)
	# 	gtable = self.model_instance.mock.galaxy_table
	# 	###TODO
		
	# 	## find the velocity of the host halo vz
	# 	## populate mock, haloID, or halo property. for each sat I need the host halo
	# 	## subtract the host halo velocity from sat velocity along oz. 
	# 	## multiply by a free parameter the velocity differente( velocity of sat in halo frame) free parameter 0.5 1.5 ?
		
	# 	# print(gtable.keys())
	# 	idx = (gtable['gal_type'] == 'centrals')
	# 	gtable['gal_type'] = idx.astype(int)
		
	# 	# range_ = gtable['gal_type'] == 0
	# 	# vdisp = self.model_instance.param_dict["vdisp"]
	# 	# gtable['vz'][range_] = (gtable['vz'][range_] - gtable['halo_vz'][range_]) * vdisp + gtable['halo_vz'][range_]  

	# 	# z_cosmo = gtable['z'].copy()

	# 	# RSD shift along OZ axis.
	# 	gtable['z'] += gtable['vz'] * self.rsd_shift
	# 	gtable['z'] %= self.Lbox
		
	# 	# import os
	# 	# np.savetxt("/global/cscratch1/sd/avariu/FastPM_SLICS/galaxy/"+os.path.basename(halo_name)+".gal.txt", np.array([gtable['x'], gtable['y'], gtable['z'], z_cosmo]).T)

	# 	sample1 = np.vstack((gtable['x'], gtable['y'], gtable['z'])).T
	# 	if self.verbose:
	# 		print('INFO: The numpy shape of the galaxy sample is {}'.format(sample1.shape))

	# 	return_dict_gal_type[halo_name + str(seed)] = gtable['gal_type']
	# 	return_dict[os.path.basename(halo_name) + str(seed)] = sample1