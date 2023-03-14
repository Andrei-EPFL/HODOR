import os
import sys
import glob
import configparser
import numpy as np
import h5py

from halotools.sim_manager import UserSuppliedHaloCatalog

class DataClass():
	""" The class that creates a list with the data to be fitted and the halo catalogs """
	def __init__(self, config_file, halo_files=None):
		config = configparser.ConfigParser()
		config.read(config_file)

		self.case_ = config['params'].getint('case')

		self.halo_mass_file_template = config['params']['halo_mass_file']
		self.halo_file_template      = config['params']['halo_file']

		self.reference_cf_template = config['params']['reference_cf']
		self.usecols_cf = tuple(map(int, config.get('params', 'usecols_cf').split(', ')))

		self.reference_pk_template = config['params']['reference_pk']
		self.usecols_pk = tuple(map(int, config.get('params', 'usecols_pk').split(', ')))

		self.reference_dens_file = config['params']['reference_dens']
		self.usecols_ndens_std = tuple(map(int, config.get('params', 'usecols_dens').split(', ')))

		self.cf_multipole = config['params'].getint('cf_multipole')
		self.pk_multipole = config['params'].getint('pk_multipole')

		if self.cf_multipole == -1 and self.pk_multipole == -1:
			print("ERROR: You have to chose at least on type of multipole to be different than -1")
			sys.exit(1)
		
		self.fit_kmin = config['pspec'].getfloat('fit_kmin')
		self.fit_kmax = config['pspec'].getfloat('fit_kmax')
		self.s_min = config['TPCF'].getfloat("s_min")
		self.s_max = config['TPCF'].getfloat("s_max")

		self.redshift  = config['params'].getfloat('redshift')
		self.box_size  = config['params'].getfloat('box_size')
		self.particle_mass = config['params'].getfloat('particle_mass')

		ndens_avg_std = np.loadtxt(self.reference_dens_file, usecols=self.usecols_ndens_std, unpack=True)
		self.ndens_avg = ndens_avg_std[0]
		print(f"INFO: The average number density of SLICS is {self.ndens_avg} +- {ndens_avg_std[1]}")

		self.LOS_list = []
		self.NLOS = 0
		if self.case_ == 0:
			if config['params']['NLOS'] == "allxslics":
				case = 21
				self.NLOS = np.inf
				self.halo_files, self.halo_mass_files, self.s_ref, self.k_ref, self.ref_clustering_vector = self.allxslics()
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
					self.halo_files, self.halo_mass_files, self.s_ref, self.k_ref, self.ref_clustering_vector = self.nlos_loslist()

				elif len(self.LOS_list) == 0:
					case = 22
					self.halo_files, self.halo_mass_files, self.s_ref, self.k_ref, self.ref_clustering_vector = self.allxslics(limit=self.NLOS)

				elif self.NLOS != len(self.LOS_list):
					print("ERROR: The NLOS is different than the size of the LOS list")
					sys.exit(1)
			print("INFO: The list of LOS is the following: ", self.LOS_list)

			self.ref_data_vector = np.append(self.ref_clustering_vector, self.ndens_avg)
		
		elif self.case_ == 1:
			self.halo_files = halo_files 
		
		self.halo_cat_list = self.compute_halocat()

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
		verify = 0
		for los in self.LOS_list:
			if os.path.isfile(self.reference_cf_template.format(los)) and os.path.isfile(self.halo_file_template.format(los)) \
			 and os.path.isfile(self.reference_pk_template.format(los)) and os.path.isfile(self.halo_mass_file_template.format(los)):

				verify = 1
				ret_halo_files.append(self.halo_file_template.format(los))
				ret_halo_mass_files.append(self.halo_mass_file_template.format(los))
				
				s, xi0, xi2, xi4 = np.loadtxt(self.reference_cf_template.format(los), usecols=self.usecols_cf, unpack=True)
				xi0_list.append(xi0)
				xi2_list.append(xi2)
				xi4_list.append(xi4)
				
				k, pk0, pk2, pk4 = np.loadtxt(self.reference_pk_template.format(los), usecols=self.usecols_pk, unpack=True)
				pk0_list.append(pk0)
				pk2_list.append(pk2)
				pk4_list.append(pk4)
				counter += 1
				print(f"INFO: Counter={counter}; LOS={los}")

		if verify == 0:
			files_ = ""
			for los in self.LOS_list:
				if not os.path.isfile(self.reference_cf_template.format(los)):
					files_ += " " + self.reference_cf_template.format(los)
				if not os.path.isfile(self.halo_file_template.format(los)):
					files_ += " " + self.halo_file_template.format(los)
				if not os.path.isfile(self.reference_pk_template.format(los)):
					files_ += " " + self.reference_pk_template.format(los)
				if not os.path.isfile(self.halo_mass_file_template.format(los)):
					files_ += " " + self.halo_mass_file_template.format(los)

			print("ERROR: No reference file or no halo or halo mass files were found. Please check the configuration file. These files do not exist: ", files_)
			sys.exit(1)

		range_s = (self.s_min <= s) & (s <= self.s_max)
		range_k = (self.fit_kmin <= k) & (k <= self.fit_kmax)

		if self.cf_multipole != -1:
			mean_xi0 = np.mean(np.array(xi0_list), 0)[range_s]
			mean_xi2 = np.mean(np.array(xi2_list), 0)[range_s]
			mean_xi4 = np.mean(np.array(xi4_list), 0)[range_s]

		if self.pk_multipole != -1:
			mean_pk0 = np.mean(np.array(pk0_list), 0)[range_k]
			mean_pk2 = np.mean(np.array(pk2_list), 0)[range_k]
			mean_pk4 = np.mean(np.array(pk4_list), 0)[range_k]

		if self.cf_multipole == -1:
			cf_data_vector = np.empty(0)
		elif self.cf_multipole == 0:
			cf_data_vector = mean_xi0
		elif self.cf_multipole == 2:
			cf_data_vector = np.concatenate((mean_xi0, mean_xi2))
		elif self.cf_multipole == 4:
			cf_data_vector = np.concatenate((mean_xi0, mean_xi2, mean_xi4))
		else:
			print("ERROR: The option for the cf_multipole is incorrect. Please use -1, 0, 2 or 4 as integers.")
			sys.exit(1)
	
		if self.pk_multipole == -1:
			pk_data_vector = np.empty(0)
		elif self.pk_multipole == 0:
			pk_data_vector = mean_pk0
		elif self.pk_multipole == 2:
			pk_data_vector = np.concatenate((mean_pk0, mean_pk2))
		elif self.pk_multipole == 4:
			pk_data_vector = np.concatenate((mean_pk0, mean_pk2, mean_pk4))
		else:
			print("ERROR: The option for the pk_multipole is incorrect. Please use -1, 0, 2 or 4 as integers.")
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
			if os.path.isfile(self.reference_cf_template.format(los)) and os.path.isfile(self.reference_pk_template.format(los)) and os.path.isfile(self.halo_mass_file_template.format(los)):
				ret_halo_files.append(file_)
				ret_halo_mass_files.append(self.halo_mass_file_template.format(los))

				s, xi0, xi2, xi4 = np.loadtxt(self.reference_cf_template.format(los), usecols=self.usecols_cf, unpack=True)
				xi0_list.append(xi0)
				xi2_list.append(xi2)
				xi4_list.append(xi4)

				k, pk0, pk2, pk4 = np.loadtxt(self.reference_pk_template.format(los), usecols=self.usecols_pk, unpack=True)
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

		if self.cf_multipole == -1:
			cf_data_vector = np.empty(0)
		elif self.cf_multipole == 0:
			cf_data_vector = mean_xi0
		elif self.cf_multipole == 2:
			cf_data_vector = np.concatenate((mean_xi0, mean_xi2))
		elif self.cf_multipole == 4:
			cf_data_vector = np.concatenate((mean_xi0, mean_xi2, mean_xi4))
		else:
			print("ERROR: The option for the cf_multipole is incorrect. Please use -1, 0, 2 or 4 as integers.")
			sys.exit(1)
	
		if self.pk_multipole == -1:
			pk_data_vector = np.empty(0)
		elif self.pk_multipole == 0:
			pk_data_vector = mean_pk0
		elif self.pk_multipole == 2:
			pk_data_vector = np.concatenate((mean_pk0, mean_pk2))
		elif self.pk_multipole == 4:
			pk_data_vector = np.concatenate((mean_pk0, mean_pk2, mean_pk4))
		else:
			print("ERROR: The option for the pk_multipole is incorrect. Please use -1, 0, 2 or 4 as integers.")
			sys.exit(1)
	
		return ret_halo_files, ret_halo_mass_files, s[range_s], k[range_k], np.concatenate((cf_data_vector, pk_data_vector))

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
			# Making all X,Y,Z values within 0-box_size
			X_data[X_data < 0] = X_data[X_data < 0] + self.box_size
			X_data[X_data > self.box_size] = X_data[X_data > self.box_size] - self.box_size
			Y_data[Y_data < 0] = Y_data[Y_data < 0] + self.box_size
			Y_data[Y_data > self.box_size] = Y_data[Y_data > self.box_size] - self.box_size
			Z_data[Z_data < 0] = Z_data[Z_data < 0] + self.box_size
			Z_data[Z_data > self.box_size] = Z_data[Z_data > self.box_size] - self.box_size
			print('INFO: min and max of X: ', np.min(X_data), np.max(X_data))
			print('INFO: min and max of Y: ', np.min(Y_data), np.max(Y_data))
			print('INFO: min and max of Z: ', np.min(Z_data), np.max(Z_data))
			# Initialize the halo catalog
			halocat = UserSuppliedHaloCatalog(redshift=self.redshift, Lbox=self.box_size, \
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