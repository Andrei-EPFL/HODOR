import os
import sys
import configparser
import numpy as np
import time
from scipy.special import erf
import warnings


def MWCens(logmass, param):
	logMcut = param[0]
	sigma_logM = param[1]
	logM1 = param[2]
	k = param[3]
	alpha = param[4]

	logM = logmass
	inv_sqrt2 = 0.7071067811865475244
	ln10 = 2.30258509299404568402
	mean_ncen = 0.5 * (1.0 + erf((ln10 * logM - ln10 * logMcut) * inv_sqrt2 / sigma_logM))
	mean_ncen[logmass < 10.75] = 0
	return mean_ncen

def MWSats(logmass, param):
	logMcut = param[0]
	sigma_logM = param[1]
	logM1 = param[2]
	k = param[3]
	alpha = param[4]

	Mcut = 10.**logMcut
	M1 = 10.**logM1
	mass = 10 ** logmass
	
	mean_nsat = np.zeros_like(mass)
	idx_nonzero = np.where(mass - k * Mcut > 0)[0]
	with warnings.catch_warnings():
	  warnings.simplefilter("ignore", RuntimeWarning)
	  mean_nsat[idx_nonzero] = ((mass[idx_nonzero] - k * Mcut) / M1)**alpha

	return mean_nsat

class DensClass():
	""" The class that contains and computes the chi2 given a model,
		the inverse of covariance matrix, the data and the config file and
		the instance of the 2pcf class """

	def __init__(self, config_file, halo_mass_files):
		config = configparser.ConfigParser()
		config.read(config_file)
		self.halo_mass_files = halo_mass_files

		self.Num_ptcl_requirement = config['params'].getfloat('Num_ptcl_requirement')
		self.particle_mass = config['params'].getfloat('particle_mass')
		self.low_lim_mass = np.log10(self.Num_ptcl_requirement * self.particle_mass)

		self.halo_mass, self.nhalos, self.stdhalos = self.compute_avg_n_halos()

	def compute_avg_n_halos(self):
		nmass_list = []
		for file_ in self.halo_mass_files:
			mass, nmass = np.loadtxt(file_, usecols=(0, 1), unpack=True)
			nmass_list.append(nmass)
		
		return mass, np.mean(np.array(nmass_list), 0), np.std(np.array(nmass_list), 0) 

	def estimate_ndens(self, param):
		mean_cens = MWCens(self.halo_mass, param)
		mean_cens[self.halo_mass < self.low_lim_mass] = 0

		mean_sats = MWSats(self.halo_mass, param)
		mean_sats[self.halo_mass < self.low_lim_mass] = 0

		return np.nansum(mean_cens * self.nhalos) + np.nansum(mean_sats * self.nhalos)