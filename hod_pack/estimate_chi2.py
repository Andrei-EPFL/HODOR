import os
import sys
import configparser
import numpy as np
import time

def compute_avg_density(dict_of_gsamples, boxL):
	counter = 0
	sum_ = 0
	for key in dict_of_gsamples.keys():
		sum_ += len(dict_of_gsamples[key])
		counter += 1
	
	avg_ = sum_ / counter
	return avg_ / (boxL ** 3)

class Chi2Class():
	""" The class that contains and computes the chi2 given a model,
		the inverse of covariance matrix, the data and the config file and
		the instance of the 2pcf class """

	def __init__(self, config_file, model_var, ref_data_vector, icov, compute_2pcf_instance, compute_pspec_instance, dens_instance):
		config = configparser.ConfigParser()
		config.read(config_file)
		self.model = config['params'].getint('model')
		self.cfmultipole = config['params'].getint('cfmultipole')
		self.pkmultipole = config['params'].getint('pkmultipole')
		self.verbose = config['params'].getboolean('verbose')
		self.outpathcf = config['params']['outpathcf']
		self.outpathcat = config['params']['outpathcat']
		if not os.path.isdir(self.outpathcf):
			print("ERROR: The output path for the CF: " + self.outpathcf + " does not exist!")
			sys.exit(1)
		self.icov = icov
		self.model_var = model_var
		self.parameters_names = model_var.parameters_names
		self.ref_data_vector = ref_data_vector
		self.compute_2pcf_instance = compute_2pcf_instance
		self.compute_pspec_instance = compute_pspec_instance
		self.dens_instance = dens_instance
		self.index_chi2 = 0

	def save_galaxy_catalogs(self, param, chunkname):
		""" Saves the galaxy catalogs given the parameters"""
		partname = chunkname + "_" + str(self.index_chi2)
		_, _, _, dict_of_gsamples, _ = self.compute_chi2_cf(param)
		for j, key in enumerate(dict_of_gsamples.keys()):
			np.savetxt(self.outpathcat  + "/CAT_" + partname + "seed" + str(key) + ".txt", dict_of_gsamples[key])
		del dict_of_gsamples

	def return_chi2_save_clustering(self, param, chunkname):
		""" Computes the chi2 and the CF given the parameters, then it saves the CF"""
		partname = chunkname + "_" + str(self.index_chi2)
		[s, xi0, xi2, xi4], [k, pk0, pk2, pk4], chisq, _, dict_gal_type = self.compute_chi2_cf(param)
		nsat, ncen = self.model_var.return_galaxy_types(dict_gal_type)
		np.savetxt(self.outpathcf + "/CF_" + partname + ".txt", np.array([s, xi0, xi2, xi4]).T)
		np.savetxt(self.outpathcf + "/PK_" + partname + ".txt", np.array([k, pk0, pk2, pk4]).T)
		del dict_gal_type
		return chisq, nsat, ncen

	def compute_best_fit(self, param):
		""" Computes the bestfit CF given the bestfit parameters"""
		[s, xi0, xi2, xi4], [k, pk0, pk2, pk4], chisq, _, _ = self.compute_chi2_cf(param)
		return [s, xi0, xi2, xi4], [k, pk0, pk2, pk4], chisq


	def compute_chi2(self, param):
		""" Computes the chi2 given the parameters"""

		_, _, chisq, _, _ = self.compute_chi2_cf(param)
		
		return chisq


	def compute_chi2_cf(self, param):
		""" Compute the chi2 given the parameters
			and returns both chi2 and the correlation function
		"""
		print('INFO: Index={} begin, param {} <<<<<<<<<<\n\n'.format(self.index_chi2 + 1, param), flush=True)

		start_time = time.time()
		params = np.atleast_1d(param)
		# if self.model == 1:
		# 	if len(params) != 9:
		# 		print("ERROR: model 1 requires 9 parameters!")
		# 		sys.exit(1)

		s, xi0, xi2, xi4 = 0, 0, 0, 0
		k, pk0, pk2, pk4 = 0, 0, 0, 0
		nsat, ncen = 0, 0
		dict_of_gsamples, dict_gal_type = 0, 0

		### Estimated using HOD curve the number density
		estimated_avg_density = self.dens_instance.estimate_ndens(param) / (self.model_var.Lbox ** 3)
		
		if np.abs(self.ref_data_vector[-1] - estimated_avg_density) < self.ref_data_vector[-1] * 0.30:
			# If the estimated number density is within 
			# The halo catalog is populated with galaxies given the parameter param.
			
			dict_of_gsamples, dict_gal_type = self.model_var.populate_mock(param,  self.ref_data_vector[-1])

			# Given the samples of galaxies, compute the 2PCF
			avg_ndens = compute_avg_density(dict_of_gsamples, self.model_var.Lbox)

			if np.abs(self.ref_data_vector[-1] - avg_ndens) < 10 * np.sqrt(1. / self.icov[-1, -1]): # Checking whether the avg density is 10 Sigma_20 larger than the reference. Sigma_20 is the error for the average of 20 catalogs.
				if self.cfmultipole == -1:
					xi = np.empty(0)
				else:
					s, xi0, xi2, xi4 = self.compute_2pcf_instance.compute_2pcf_corrfunc(dict_of_gsamples, self.index_chi2)
					
					if self.cfmultipole == 0:
						xi = xi0
					elif self.cfmultipole == 2:
						xi = np.concatenate([xi0, xi2])
					elif self.cfmultipole == 4:
						xi = np.concatenate([xi0, xi2, xi4])

				if self.pkmultipole == -1:
					pk = np.empty(0)
				else:
					k, pk0, pk2, pk4 = self.compute_pspec_instance.pspec_2(dict_of_gsamples, self.index_chi2)	
					if self.pkmultipole == 0:
						pk = pk0
					elif self.pkmultipole == 2:
						pk = np.concatenate([pk0, pk2])
					elif self.pkmultipole == 4:
						pk = np.concatenate([pk0, pk2, pk4])

				model_vector = np.hstack((xi, pk, avg_ndens))
				print(model_vector.shape)
				if self.verbose:
					print("STATUS: Computing chi2...")

				diff = model_vector - self.ref_data_vector
				matrix = self.icov

				chisq = np.dot(diff.T, np.dot(matrix, diff))
			else:
				chisq = 10e9 * self.icov[-1, -1] * (avg_ndens - self.ref_data_vector[-1])**2 
				# chisq = 1e101 
		
			nsat, ncen = self.model_var.return_galaxy_types(dict_gal_type)

		else:
			chisq = 10e9 * self.icov[-1, -1] * (estimated_avg_density - self.ref_data_vector[-1])**2 
			# chisq = 1e101
		
		self.index_chi2 = self.index_chi2 + 1
		end_time = time.time()

		print('INFO: Index={} end, chi2 = {} in {} seconds ncen {} nsat {} param {} <<<<<<<<<<\n\n'.format(self.index_chi2, chisq, end_time - start_time, ncen, nsat, param), flush=True)

		return [s, xi0, xi2, xi4], [k, pk0, pk2, pk4], chisq, dict_of_gsamples, dict_gal_type
