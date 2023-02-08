import os
import sys
import configparser
import numpy as np
import time

def compute_avg_density(dict_of_gsamples, box_size):
	counter = 0
	sum_ = 0
	for key in dict_of_gsamples.keys():
		sum_ += len(dict_of_gsamples[key])
		counter += 1
	
	avg_ = sum_ / counter
	return avg_ / (box_size ** 3)

class Chi2Class():
	""" The class that contains and computes the chi2 given a model,
		the inverse of covariance matrix, the data and the config file and
		the instance of the 2pcf class """

	def __init__(self, config_file, model_instance, ref_data_vector, icov, compute_2pcf_instance, compute_pspec_instance, dens_instance):
		config = configparser.ConfigParser()
		config.read(config_file)

		self.model        = config['params'].getint('model')
		self.cf_multipole = config['params'].getint('cf_multipole')
		self.pk_multipole = config['params'].getint('pk_multipole')
		self.verbose      = config['params'].getboolean('verbose')

		self.outpath_clustering = config['params']['outpath_clustering']
		self.outpath_catalog    = config['params']['outpath_catalog']

		if not os.path.isdir(self.outpath_clustering):
			print("ERROR: The output path for the CF: " + self.outpath_clustering + " does not exist!")
			sys.exit(1)

		self.icov = icov
		self.model_instance   = model_instance
		self.parameters_names = model_instance.parameters_names
		self.ref_data_vector  = ref_data_vector

		self.compute_2pcf_instance  = compute_2pcf_instance
		self.compute_pspec_instance = compute_pspec_instance
		self.dens_instance          = dens_instance

		self.index_chi2 = 0


	def save_galaxy_catalogs(self, param, chunkname):
		""" Saves the galaxy catalogs given the parameters"""
		partname = chunkname + "_" + str(self.index_chi2)
		_, _, _, dict_of_gsamples = self.compute_chi2_cf(param)
		for j, key in enumerate(dict_of_gsamples.keys()):
			np.savetxt(f"{self.outpath_catalog}/CAT_{partname}seed{str(key)}.txt", dict_of_gsamples[key])
		del dict_of_gsamples


	def return_chi2_save_clustering(self, param, chunkname):
		""" Computes the chi2 and the CF given the parameters, then it saves the CF"""
		partname = chunkname + "_" + str(self.index_chi2)
		[s, xi0, xi2, xi4], [k, pk0, pk2, pk4], chi2, dict_of_gsamples = self.compute_chi2_cf(param)
		nsat, ncen = self.model_instance.return_galaxy_types(dict_of_gsamples)
		np.savetxt(f"{self.outpath_clustering}/CF_{partname}.txt", np.array([s, xi0, xi2, xi4]).T)
		np.savetxt(f"{self.outpath_clustering}/PK_{partname}.txt", np.array([k, pk0, pk2, pk4]).T)
		return chi2, nsat, ncen


	def compute_best_fit(self, param):
		""" Computes the bestfit CF given the bestfit parameters"""
		[s, xi0, xi2, xi4], [k, pk0, pk2, pk4], chi2, dict_of_gsamples = self.compute_chi2_cf(param)
		nsat, ncen = self.model_instance.return_galaxy_types(dict_of_gsamples)
		return [s, xi0, xi2, xi4], [k, pk0, pk2, pk4], chi2, [nsat, ncen]


	def compute_chi2(self, param):
		""" Computes the chi2 given the parameters"""
		_, _, chi2, _ = self.compute_chi2_cf(param)
		return chi2


	def compute_chi2_cf(self, param):
		""" Compute the chi2 given the parameters
			and returns both chi2 and the correlation function
		"""
		print('INFO: Index={} begin, param {} <<<<<<<<<<\n\n'.format(self.index_chi2, param), flush=True)

		start_time = time.time()
		params = np.atleast_1d(param)

		### Initial values
		s, xi0, xi2, xi4 = 0, 0, 0, 0
		k, pk0, pk2, pk4 = 0, 0, 0, 0
		nsat, ncen = 0, 0
		dict_of_gsamples = 0

		### Estimated the number density using the HOD curve 
		estimated_avg_density = self.dens_instance.estimate_ndens(param) / (self.model_instance.box_size ** 3)
		
		if np.abs(self.ref_data_vector[-1] - estimated_avg_density) < self.ref_data_vector[-1] * 0.30:
			# If the estimated number density is within 
			# The halo catalog is populated with galaxies given the parameter param.
			
			dict_of_gsamples = self.model_instance.populate_mock(param, self.ref_data_vector[-1])

			# Given the samples of galaxies, compute the average density
			avg_ndens = compute_avg_density(dict_of_gsamples, self.model_instance.box_size)

			if np.abs(self.ref_data_vector[-1] - avg_ndens) < 10 * np.sqrt(1. / self.icov[-1, -1]): # Checking whether the avg density is 10 Sigma_20 larger than the reference. Sigma_20 is the error for the average of 20 catalogs.
				if self.cf_multipole == -1:
					xi = np.empty(0)
				else:
					s, xi0, xi2, xi4 = self.compute_2pcf_instance.compute_2pcf_pyfcfc(dict_of_gsamples, self.index_chi2)

					if self.cf_multipole == 0:
						xi = xi0
					elif self.cf_multipole == 2:
						xi = np.concatenate([xi0, xi2])
					elif self.cf_multipole == 4:
						xi = np.concatenate([xi0, xi2, xi4])

				if self.pk_multipole == -1:
					pk = np.empty(0)
				else:
					k, pk0, pk2, pk4 = self.compute_pspec_instance.pspec(dict_of_gsamples, self.index_chi2)	

					if self.pk_multipole == 0:
						pk = pk0
					elif self.pk_multipole == 2:
						pk = np.concatenate([pk0, pk2])
					elif self.pk_multipole == 4:
						pk = np.concatenate([pk0, pk2, pk4])

				model_vector = np.hstack((xi, pk, avg_ndens))

				print(f"INFO: The shape of the model vector is {model_vector.shape}")
				if self.verbose:
					print("STATUS: Computing chi2...")

				diff = model_vector - self.ref_data_vector
				matrix = self.icov

				chi2 = np.dot(diff.T, np.dot(matrix, diff))
			else:
				chi2 = 10e9 * self.icov[-1, -1] * (avg_ndens - self.ref_data_vector[-1])**2 
				# chi2 = 1e101 
		
			nsat, ncen = self.model_instance.return_galaxy_types(dict_of_gsamples)

		else:
			chi2 = 10e9 * self.icov[-1, -1] * (estimated_avg_density - self.ref_data_vector[-1])**2 
			# chi2 = 1e101

		end_time = time.time()

		print('INFO: Index={} end, chi2 = {} in {} seconds ncen {} nsat {} param {} <<<<<<<<<<\n\n'.format(self.index_chi2, chi2, end_time - start_time, ncen, nsat, param), flush=True)

		self.index_chi2 = self.index_chi2 + 1
		return [s, xi0, xi2, xi4], [k, pk0, pk2, pk4], chi2, dict_of_gsamples
