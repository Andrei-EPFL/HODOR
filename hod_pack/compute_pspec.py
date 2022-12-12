import configparser
from ctypes import *
import numpy as np


class Data(Structure):
	_fields_ = [
		("x", c_double * 3),
		("w", c_double)
	]


class Cata(Structure):
	_fields_ = [
		("num", c_int), # number of data catalogs to be read
		("data", POINTER(POINTER(Data))),
		("rand", POINTER(POINTER(Data))),
		("ndata", POINTER(c_size_t)),    
		("nrand", POINTER(c_size_t)),
		("wdata", POINTER(c_double)),
		("wrand", POINTER(c_double)),
		("alpha", POINTER(c_double)),
		("shot", POINTER(c_double)),
		("norm", POINTER(c_double)),
	]


class ComputePKClass():
	""" The class that computes the pspec given the config file and a sample of galaxies """
	def __init__(self, config_file):
		config = configparser.ConfigParser()
		config.read(config_file)

		self.fit_kmin = config['pspec'].getfloat('fit_kmin')
		self.fit_kmax = config['pspec'].getfloat('fit_kmax')
		self.Lbox = config['params'].getfloat('box_size')


	def compute_catalog(self, x, y, z):

		x = x % self.Lbox
		y = y % self.Lbox
		z = z % self.Lbox

		num = 1
		ndata = len(x)
		nrand = 0
		wdata = ndata
		wrand = 0.
		alpha = 0.
		shot = 0.
		norm = 0.


		data_instance = Data * ndata
		data_array = data_instance()

		cata_instance = Cata
		ca = cata_instance()

		ca.num = c_int(num)
		
		ca.data.contents = pointer(data_array)
		ca.rand= None
		
		ca.ndata.contents = c_size_t(ndata)
		ca.nrand.contents = c_double(nrand)

		ca.wdata.contents = c_double(wdata)
		ca.wrand.contents = c_double(wrand)
		
		ca.alpha.contents = c_double(alpha)
		ca.shot.contents = c_double(shot)
		ca.norm.contents = c_double(norm)

		for i in range(ndata):
			da = ca.data[0][i]

			da.w = c_double(1.)
			da.x[0] = c_double(x[i])
			da.x[1] = c_double(y[i])
			da.x[2] = c_double(z[i])
	
	
		return ca


	def pspec(self, dict_of_gsamples, index):
		lib = cdll.LoadLibrary("/global/homes/a/avariu/phd/chengscodes/powspec_cffi/libpowspec.so")
		config_file = b'/global/homes/a/avariu/phd/chengscodes/powspec_cffi/etc/powspec_HODFIT.conf'
		output_file = b'./temp_to_delete.temp'
		input_file =  b'/global/homes/a/avariu/phd/chengscodes/powspec_cffi/etc/powspec_HODFIT.conf'

		### Arguments
		size_text = len(config_file)
		arg0  = create_string_buffer(b'test', size_text)
		arg1  = create_string_buffer(b'-c', size_text)
		arg2 = create_string_buffer(config_file)

		arg3  = create_string_buffer(b'-d', size_text)
		arg4  = create_string_buffer(input_file, size_text)

		arg5  = create_string_buffer(b'-a', size_text)
		arg6  = create_string_buffer(output_file, size_text)

		N_args = 7

		args_pointer_instance = POINTER(c_char) * N_args
		args_pointer_arr_to_send = args_pointer_instance(arg0, arg1, arg2, arg3, arg4, arg5, arg6)

		### Multipoles list
		pk0_list = []
		pk2_list = []
		pk4_list = []
		
		### Go through all galaxy samples
		for j, key in enumerate(dict_of_gsamples.keys()):
			cata = self.compute_catalog(dict_of_gsamples[key]["x"], dict_of_gsamples[key]["y"], dict_of_gsamples[key]["z_rsd"])

			nkbin = c_int(0)
			nbin = 1000
			pk_instance = c_double * (4 * nbin)
			pk_array = pk_instance()

			lib.compute_pk.restype = c_int

			value_return = lib.compute_pk(pointer(cata), pointer(nkbin), pointer(pk_array), c_int(N_args), args_pointer_arr_to_send)
			if value_return == 0:
				print("ERROR: the pk code crashed")
				del pk_array
				del cata
				return 0, 0, 0, 0

			nbin = nkbin.value

			### From C to NumPy Array
			k = np.zeros(nbin)
			pk0 = np.zeros(nbin)
			pk2 = np.zeros(nbin)
			pk4 = np.zeros(nbin)
			
			for i in range(nbin):
				k[i]   = pk_array[i]
				pk0[i] = pk_array[nbin + i]
				pk2[i] = pk_array[2 * nbin + i]
				pk4[i] = pk_array[3 * nbin + i]

			del pk_array
			# del cata.data
			del cata
			## import gc
			## gc.collect()
			pk0_list.append(pk0)
			pk2_list.append(pk2)
			pk4_list.append(pk4)

		range_ = (self.fit_kmin <= k) & (k <= self.fit_kmax)

		mean_pk0 = np.mean(np.array(pk0_list), 0)[range_]
		mean_pk2 = np.mean(np.array(pk2_list), 0)[range_]
		mean_pk4 = np.mean(np.array(pk4_list), 0)[range_]

		return k[range_], mean_pk0, mean_pk2, mean_pk4

