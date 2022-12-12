import configparser
import numpy as np
from halotools.mock_observables import s_mu_tpcf
from halotools.mock_observables import tpcf_multipole
import Corrfunc.theory as cft


class ComputeCFClass():
	""" The class that computes the 2PCF given the config file and a sample of galaxies """
	def __init__(self, config_file):
		config = configparser.ConfigParser()
		config.read(config_file)

		self.verbose = config['params'].getboolean('verbose')
		self.s_min = config['TPCF'].getfloat("s_min")
		self.s_max = config['TPCF'].getfloat("s_max")
		self.n_s_bins = config['TPCF'].getint("n_s_bins")
		#self.s_bins = np.geomspace(self.s_min, self.s_max, self.n_s_bins + 1)
		self.s_bins = np.linspace(self.s_min, self.s_max, self.n_s_bins + 1)
		print("INFO: The s_bins are: ", self.s_bins)
		self.r = np.array([(a + b) / 2.0 for a, b in zip(self.s_bins[:-1], self.s_bins[1:])])

		self.mu_min = config['TPCF'].getfloat("mu_min")
		self.mu_max = config['TPCF'].getfloat("mu_max")
		self.n_mu_bins = config['TPCF'].getint("n_mu_bins")
		self.mu_bins = np.linspace(self.mu_min, self.mu_max, self.n_mu_bins + 1)
		print("INFO: The mu_bins are: ", self.mu_bins)

		self.nthreads = config['TPCF'].getint("nthreads")
		self.Lbox = config['params'].getfloat('Lbox')

	def compute_2pcf_halotools(self, dict_of_gsamples):
		""" Compute the 2PCF of galaxies given the
		dictionary of galaxies. It uses the halotools package to compute the 2PCF """

		if self.verbose:
			print('\nSTATUS: Compute the 2pcf...')

		xi0_arr = np.zeros((len(dict_of_gsamples.keys()), len(self.r)))
		xi2_arr = np.zeros((len(dict_of_gsamples.keys()), len(self.r)))
		xi4_arr = np.zeros((len(dict_of_gsamples.keys()), len(self.r)))

		for j, key in enumerate(dict_of_gsamples.keys()):
			if self.verbose:
				print('INFO: >> Seed = {}...'.format(key))
			sample1 = dict_of_gsamples[key]

			# Compute xi0, xi2, xi4
			if self.verbose:
				print('SubSTATUS: Compute the 2D 2pcf...')
			xi_s_mu = s_mu_tpcf(sample1, self.s_bins, self.mu_bins, period=self.Lbox, num_threads=self.nthreads)
			if self.verbose:
				print('SubSTATUS: Compute the multipoles...')
			xi0 = tpcf_multipole(xi_s_mu, self.mu_bins, order=0)
			xi2 = tpcf_multipole(xi_s_mu, self.mu_bins, order=2)
			xi4 = tpcf_multipole(xi_s_mu, self.mu_bins, order=4)

			xi0_arr[j] = xi0
			xi2_arr[j] = xi2
			xi4_arr[j] = xi4

		mean_xi0 = np.mean(xi0_arr, 0)
		mean_xi2 = np.mean(xi2_arr, 0)
		mean_xi4 = np.mean(xi4_arr, 0)
		return self.r, mean_xi0, mean_xi2, mean_xi4

	def spherical_sector_volume(self, s, mu):
		"""
		This function is used to calculate analytical randoms.
		Calculate the volume of a spherical sector, used for the analytical randoms.
		https://en.wikipedia.org/wiki/Spherical_sector
		Note that the extra factor of 2 is to get the reflection.
		"""
		theta = np.arccos(mu)

		vol = (2.0*np.pi/3.0) * np.outer((s**3.0), (1.0-np.cos(theta)))*2.0
		return vol

	def compute_analytical_RR(self, NR):
		""" This function is used to compute an analytical value of the
		RR term of the 2PCF """

		# do volume calculations
		mu_bins_reverse_sorted = np.sort(self.mu_bins)[::-1]
		dv = self.spherical_sector_volume(self.s_bins, mu_bins_reverse_sorted)
		dv = np.diff(dv, axis=1)  # volume of wedges
		dv = np.diff(dv, axis=0)  # volume of wedge 'pieces'
		global_volume = self.Lbox**3

		# calculate the random-random pairs.
		rhor = NR**2/global_volume
		RR = (dv*rhor)
		return RR

	def compute_2pcf_corrfunc(self, dict_of_gsamples, index):
		""" Compute the 2PCF of galaxies given the
		dictionary of galaxies. It uses the corrfunc package (which is much faster than halotools)
		to compute the 2PCF, with a natural estimator """

		if self.verbose:
			print('\nSTATUS: Compute the 2pcf...')

		xi0_arr = np.zeros((len(dict_of_gsamples.keys()), len(self.r)))
		xi2_arr = np.zeros((len(dict_of_gsamples.keys()), len(self.r)))
		xi4_arr = np.zeros((len(dict_of_gsamples.keys()), len(self.r)))

		for j, key in enumerate(dict_of_gsamples.keys()):
			if self.verbose:
				print('INFO: >> Seed = {}...'.format(key))
			sample1 = dict_of_gsamples[key]
			# Compute xi0, xi2, xi4
			if self.verbose:
				print('SubSTATUS: Compute the RR term...')
			RR = self.compute_analytical_RR(len(sample1))

			autocorr = 1
			if self.verbose:
				print('SubSTATUS: Compute the DD term...')
			DD_s_mu = cft.DDsmu(autocorr, self.nthreads, self.s_bins, self.mu_max, self.n_mu_bins, sample1[:, 0], sample1[:, 1], sample1[:, 2], boxsize=self.Lbox, verbose=self.verbose)
			DD_s_mu = np.reshape(DD_s_mu[:]['npairs'], (self.n_s_bins, self.n_mu_bins))
			# Corrfunc adds also the central galaxy in an s bin (when smin==0). So one needs to subtract it from the number of pairs.
			if self.s_min == 0:
				DD_s_mu[0][0] = DD_s_mu[0][0] - len(sample1)

			if self.verbose:
				print('SubSTATUS: Compute the 2D 2pcf...')
			xi_s_mu = (DD_s_mu / RR) - 1

			if self.verbose:
				print('SubSTATUS: Compute the multipoles...')
			xi0 = tpcf_multipole(xi_s_mu, self.mu_bins, order=0)
			xi2 = tpcf_multipole(xi_s_mu, self.mu_bins, order=2)
			xi4 = tpcf_multipole(xi_s_mu, self.mu_bins, order=4)

			xi0_arr[j] = xi0
			xi2_arr[j] = xi2
			xi4_arr[j] = xi4
			#np.savetxt("/global/homes/a/avariu/phd/hodsimulations/output/fit/scipy/cf_gal/2nd_1.087e+00_2.203e-03_-5.980e-01_1.513e+02_9.998e-01_1seed_Ginevra/CF_" + str(j) + "_" + str(key) +".txt", np.array([self.r, xi0, xi2, xi4]).T)

		mean_xi0 = np.mean(xi0_arr, 0)
		mean_xi2 = np.mean(xi2_arr, 0)
		mean_xi4 = np.mean(xi4_arr, 0)
		return self.r, mean_xi0, mean_xi2, mean_xi4
	