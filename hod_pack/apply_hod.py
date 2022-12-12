import os
import sys
import configparser
from copy import copy
import numpy as np
from halotools.empirical_models import HodModelFactory
from halotools.empirical_models import TrivialPhaseSpace, NFWPhaseSpace

from .hod_models import MWCens, MWSats, ContrerasCens, ContrerasSats

class ModelClass():
	""" The class that populates the halo catalog with galaxies
	given the config file """
	def __init__(self, config_file, halo_files, halo_cat_list):
		config = configparser.ConfigParser()
		config.read(config_file)

		self.model = config['params'].getint('model')
		if self.model == 0:
			model_name = 'MW'
		elif self.model == 1:
			model_name = 'GP18'
		else:
			raise ValueError('Unrecognised model')
		self.parameters_names = tuple(map(str, config.get('params', 'model_params_names').split(', ')))
		self.num_params = config['params'].getint('num_params')

		print(f'\nINFO: The used HOD model is {model_name} having {self.parameters_names} as parameters')
		if len(self.parameters_names) != self.num_params:
			print("ERROR: Check the num_params and model_params_names in the config file")
			sys.exit(1)

		self.verbose   = config['params'].getboolean('verbose')
		self.redshift  = config['params'].getfloat('redshift')
		self.box_size  = config['params'].getfloat('box_size')
		self.Omega_m   = config['params'].getfloat('Omega_m')
		self.init_seed = config['params'].getint('init_seed')
		self.num_seeds = config['params'].getint('num_seeds')

		print(f'\nINFO: The number of seeds to populate galaxies is {self.num_seeds}, starting from {self.init_seed}')
		self.z_space = config['params'].getboolean('z_space')
		self.halo_files = halo_files
		print(f'\nINFO: The number of halo catalogs used to populate galaxies is {len(halo_files)}')

		self.Num_ptcl_requirement = config['params'].getfloat('Num_ptcl_requirement')

		self.halo_cat_list = halo_cat_list
		self.model_instance = self.compute_model_instance()
		self.rsd_shift = self.compute_rsd_shift()


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


	def compute_rsd_shift(self):
		""" Compute the RSD redshift """
		if self.z_space:
			Omega_l = 1.0 - self.Omega_m
			hubble = 100 * np.sqrt(self.Omega_m * (self.redshift + 1.0)**3 + Omega_l)
			rsd_shift = (self.redshift + 1.0) / hubble
		else:
			rsd_shift = 0
		return rsd_shift


	def populate_mock(self, param_vals, ref_num_dens, indx=0):
		""" Populate the halo catalog with galaxies for different seeds """
		if self.verbose:
			print('\nSTATUS: Populate the halo catalog with galaxies...')
		for i in range(self.num_params):
			self.model_instance.param_dict[self.parameters_names[i]] = param_vals[i]

		### Function which tells you how many central and satellite galaxies were alocated to the FastPM halo catalog.
		#self.check_galaxy_types(self.init_seed)
		#exit()

		### Actually populating the halo mocks with galaxies for different seeds.
		return_dict = {}
		for j, halo_cat in enumerate(self.halo_cat_list):
			for i in range(self.num_seeds):
				seedt, gal_samplet = self.intermediate_populate_mock(self.init_seed + i*1000 + j + indx, halo_cat)

				return_dict[os.path.basename(self.halo_files[j]) + str(seedt)] = gal_samplet

				if (i == 0) and (j == 0):
					meas_num_dens = len(gal_samplet) / (self.box_size ** 3)
					if np.abs(meas_num_dens - ref_num_dens) > 0.20 * ref_num_dens: ### Check whether one catalog is 8 Sigma_1 larger than the reference. Sigma_1 is the error of one catalog. 
						return return_dict

		return return_dict


	def intermediate_populate_mock(self, seed, halo_cat):
		""" This is an intermediate step in populating the halo catalog with galaxies,
		required to parallelize the procedure for multiple seeds """
		if self.verbose:
			print('SubSTATUS: Populate the halo catalog with galaxies using {} as seed...'.format(seed))
		
		self.model_instance.populate_mock(halo_cat, Num_ptcl_requirement=self.Num_ptcl_requirement, seed=seed)
		gtable = self.model_instance.mock.galaxy_table
		idx = (gtable['gal_type'] == 'centrals')
		gtable['gal_type'] = idx.astype(int)

		### TO DO for velocity dispersion:
		### Must include an if, for the case when the vdisp is not used.
		### Velocity Dispersion
		range_ = gtable['gal_type'] == 0
		vdisp = self.model_instance.param_dict["vdisp"]
		gtable['vz'][range_] = (gtable['vz'][range_] - gtable['halo_vz'][range_]) * vdisp + gtable['halo_vz'][range_]

		# RSD shift along OZ axis.
		z_rsd = gtable['z'] + gtable['vz'] * self.rsd_shift
		z_rsd = (z_rsd + self.box_size) % self.box_size

		gal_sample = np.rec.fromarrays([gtable['x'], gtable['y'], z_rsd, gtable['z'], gtable['gal_type']], dtype=[('x', 'f4'), ('y', 'f4'), ('z_rsd', 'f4'), ('z', 'f4'), ('gal_type', 'i4')])

		return seed, gal_sample


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


	def return_galaxy_types(self, dict_of_gsamples):
		""" Function which tells you how many central and satellite galaxies were alocated to the FastPM halo catalog. """

		nsat = np.zeros(len(dict_of_gsamples.keys()))
		ncen = np.zeros(len(dict_of_gsamples.keys()))

		for j, key in enumerate(dict_of_gsamples.keys()):
			if self.verbose:
				print('INFO: >> Seed = {}...'.format(key))

			gal_type = dict_of_gsamples[key]["gal_type"]

			unique, frequency = np.unique(gal_type, return_counts=True)
			# print unique values array
			#print("Unique Values:", unique)
			nsat[j] = frequency[0]	

			if len(frequency) == 2:
				ncen[j] = frequency[1]
			# print frequency array
			#print("Frequency Values:", frequency)

		return np.mean(nsat), np.mean(ncen)
