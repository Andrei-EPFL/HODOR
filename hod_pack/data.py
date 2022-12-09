import os
import sys
import glob
import configparser
import numpy as np

class DataClass():
	""" The class that creates a list with the data to be fitted and the halo catalogs """
	def __init__(self, config_file):
		config = configparser.ConfigParser()
		config.read(config_file)

		self.halo_mass_file_template = config['params']['hmassfile']
		self.halo_file_template = config['params']['halofile']
		self.fit_file_template = config['params']['fitfile']
		self.pk_file_template = config['params']['pkfile']
		self.dens_fit_file = config['params']['densfile']
		
		self.cfmultipole = config['params'].getint('cfmultipole')
		self.pkmultipole = config['params'].getint('pkmultipole')

		if self.cfmultipole == -1 and self.pkmultipole == -1:
			print("ERROR: You have to chose at least on type of multipole to be different than -1")
			sys.exit(1)
		
		self.fit_kmin = config['pspec'].getfloat('fit_kmin')
		self.fit_kmax = config['pspec'].getfloat('fit_kmax')
		self.s_min = config['TPCF'].getfloat("s_min")
		self.s_max = config['TPCF'].getfloat("s_max")
		
		ndens_std = np.loadtxt(self.dens_fit_file, usecols=(0), unpack=True)
		self.ndens_avg = ndens_std[0]
		print(f"INFO: The average number density of SLICS is {self.ndens_avg} +- {ndens_std[1]}")

		self.LOS_list = []
		self.NLOS = 0

		if config['params']['NLOS'] == "allxslics":
			case = 21
			self.NLOS = np.inf
			self.halo_files, self.halo_mass_files, self.s_ref, self.k_ref, self.ref_data_vector = self.allxslics()
		else:
			try:
				self.NLOS = config['params'].getint('NLOS')
			except:
				print("ERROR: NLOS is not an int or allxslics")
				sys.exit(1)

			# try:
			self.LOS_list = tuple(map(int, config.get('params', 'LOS').split(', ')))
			# except:
				# print("Warning: the LOS parameter in config file is empty.", flush=True)
			
			if self.NLOS == len(self.LOS_list):
				case = 23
				self.halo_files, self.halo_mass_files, self.s_ref, self.k_ref, self.ref_data_vector = self.nlos_loslist()

			elif len(self.LOS_list) == 0:
				case = 22
				self.halo_files, self.halo_mass_files, self.s_ref, self.k_ref, self.ref_data_vector = self.allxslics(limit=self.NLOS)

			elif self.NLOS != len(self.LOS_list):
				print("ERROR: The NLOS is different than the size of the LOS list")
				sys.exit(1)

		print("TESTSETST", self.LOS_list)
			
		# Case 1:
		## One LOS multiple seeds

		# Case 2:
		## Multiple LOS one seed
			### Case 21)
			#### all possible LOS

			### Case 22)
			#### given number of LOS

			### Case 23)
			#### given specific LOS

		# Case 3:
		## Multiple LOS multiple seeds

	def nlos_loslist(self):
		ret_halo_files = []
		ret_halo_mass_files = []

		xi0_list = []
		xi2_list = []
		xi4_list = []
		
		pk0_list = []
		pk2_list = []
		pk4_list = []
		
		counter = 0
		
		for los in self.LOS_list:
			if os.path.isfile(self.fit_file_template.format(los)) and os.path.isfile(self.halo_file_template.format(los)) and os.path.isfile(self.pk_file_template.format(los)) and os.path.isfile(self.halo_mass_file_template.format(los)):
				ret_halo_files.append(self.halo_file_template.format(los))
				ret_halo_mass_files.append(self.halo_mass_file_template.format(los))
				
				s, xi0, xi2, xi4 = np.loadtxt(self.fit_file_template.format(los), usecols=(0,1,2,3), unpack=True)
				xi0_list.append(xi0)
				xi2_list.append(xi2)
				xi4_list.append(xi4)
				
				k, pk0, pk2, pk4 = np.loadtxt(self.pk_file_template.format(los), usecols=(0,1,2,3), unpack=True)
				pk0_list.append(pk0)
				pk2_list.append(pk2)
				pk4_list.append(pk4)
				
				counter += 1
				print(f"INFO: Counter={counter}; LOS={los}")

		range_s = (self.s_min <= s) & (s <= self.s_max)
		mean_xi0 = np.mean(np.array(xi0_list), 0)[range_s]
		mean_xi2 = np.mean(np.array(xi2_list), 0)[range_s]
		mean_xi4 = np.mean(np.array(xi4_list), 0)[range_s]

		range_k = (self.fit_kmin <= k) & (k <= self.fit_kmax)
		mean_pk0 = np.mean(np.array(pk0_list), 0)[range_k]
		mean_pk2 = np.mean(np.array(pk2_list), 0)[range_k]
		mean_pk4 = np.mean(np.array(pk4_list), 0)[range_k]

		if self.cfmultipole == -1:
			cf_data_vector = np.empty(0)
		elif self.cfmultipole == 0:
			cf_data_vector = mean_xi0
		elif self.cfmultipole == 2:
			cf_data_vector = np.concatenate((mean_xi0, mean_xi2))
		elif self.cfmultipole == 4:
			cf_data_vector = np.concatenate((mean_xi0, mean_xi2, mean_xi4))
		else:
			print("ERROR: The option for the cfmultipole is incorrect. Please use -1, 0, 2 or 4 as integers.")
			sys.exit(1)
	
		if self.pkmultipole == -1:
			pk_data_vector = np.empty(0)
		elif self.pkmultipole == 0:
			pk_data_vector = mean_pk0
		elif self.pkmultipole == 2:
			pk_data_vector = np.concatenate((mean_pk0, mean_pk2))
		elif self.pkmultipole == 4:
			pk_data_vector = np.concatenate((mean_pk0, mean_pk2, mean_pk4))
		else:
			print("ERROR: The option for the pkmultipole is incorrect. Please use -1, 0, 2 or 4 as integers.")
			sys.exit(1)
	
		return ret_halo_files, ret_halo_mass_files, s[range_s], k[range_k], np.concatenate((cf_data_vector, pk_data_vector))

	def allxslics(self, limit=np.inf):
		ret_halo_files = []
		ret_halo_mass_files = []

		xi0_list = []
		xi2_list = []
		xi4_list = []

		pk0_list = []
		pk2_list = []
		pk4_list = []
		
		counter = 0
		halo_files = glob.glob(self.halo_file_template.format("*"))
		print(f"INFO: There are {len(halo_files)} halo files.")
		for file_ in halo_files:
			beg = file_.find("LOS")
			end = file_.find("fof")
			los = file_[beg + 3: end - 1]
			if os.path.isfile(self.fit_file_template.format(los)) and os.path.isfile(self.pk_file_template.format(los)) and os.path.isfile(self.halo_mass_file_template.format(los)):
				ret_halo_files.append(file_)
				ret_halo_mass_files.append(self.halo_mass_file_template.format(los))

				s, xi0, xi2, xi4 = np.loadtxt(self.fit_file_template.format(los), usecols=(0,1,2,3), unpack=True)
				xi0_list.append(xi0)
				xi2_list.append(xi2)
				xi4_list.append(xi4)

				k, pk0, pk2, pk4 = np.loadtxt(self.pk_file_template.format(los), usecols=(0,1,2,3), unpack=True)
				pk0_list.append(pk0)
				pk2_list.append(pk2)
				pk4_list.append(pk4)

				counter += 1
				print(f"INFO: Counter={counter}; LOS={los}")

			if counter >= limit:
				break
				
		range_s = (self.s_min <= s) & (s <= self.s_max)
		mean_xi0 = np.mean(np.array(xi0_list), 0)[range_s]
		mean_xi2 = np.mean(np.array(xi2_list), 0)[range_s]
		mean_xi4 = np.mean(np.array(xi4_list), 0)[range_s]

		range_k = (self.fit_kmin <= k) & (k <= self.fit_kmax)
		mean_pk0 = np.mean(np.array(pk0_list), 0)[range_k]
		mean_pk2 = np.mean(np.array(pk2_list), 0)[range_k]
		mean_pk4 = np.mean(np.array(pk4_list), 0)[range_k]

		if self.cfmultipole == -1:
			cf_data_vector = np.empty(0)
		elif self.cfmultipole == 0:
			cf_data_vector = mean_xi0
		elif self.cfmultipole == 2:
			cf_data_vector = np.concatenate((mean_xi0, mean_xi2))
		elif self.cfmultipole == 4:
			cf_data_vector = np.concatenate((mean_xi0, mean_xi2, mean_xi4))
		else:
			print("ERROR: The option for the cfmultipole is incorrect. Please use -1, 0, 2 or 4 as integers.")
			sys.exit(1)
	
		if self.pkmultipole == -1:
			pk_data_vector = np.empty(0)
		elif self.pkmultipole == 0:
			pk_data_vector = mean_pk0
		elif self.pkmultipole == 2:
			pk_data_vector = np.concatenate((mean_pk0, mean_pk2))
		elif self.pkmultipole == 4:
			pk_data_vector = np.concatenate((mean_pk0, mean_pk2, mean_pk4))
		else:
			print("ERROR: The option for the pkmultipole is incorrect. Please use -1, 0, 2 or 4 as integers.")
			sys.exit(1)
	
		return ret_halo_files, ret_halo_mass_files, s[range_s], k[range_k], np.concatenate((cf_data_vector, pk_data_vector))
