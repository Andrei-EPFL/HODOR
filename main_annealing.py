""" Apply HOD and fit the 2PCF to High resolution simulation """
#!/global/homes/a/avariu/.conda/envs/hodfit/bin/python
import sys
import os
import shutil
import time
import argparse
import configparser
import numpy as np
from chi2 import Chi2Class
from apply_hod import ModelClass
from minimizer import ClassMinimizer
from compute_2pcf import ComputeCFClass
# from compute_pspec import ComputePKClass

from data import DataClass

def main():
	""" Main function of the code"""
	# Parameters
	start = time.time()
	parser = argparse.ArgumentParser()
	parser.add_argument("--config", type=str, help="ini file holding configuration")

	args = parser.parse_args()

	#### Configuration file
	config_file = args.config # config file
	if config_file is None:
		config_file = '/global/homes/a/avariu/phd/hodsimulations/HOD_pipeline_Cheng/fit/config.ini'

	if not os.path.isfile(config_file):
		print("ERROR: The config file: " + config_file + " does not exist!")
		sys.exit(1)

	config_name = os.path.basename(config_file)

	config = configparser.ConfigParser()
	config.read(config_file)

	outmultinest = config['params']['outmultinest']
	if not os.path.isdir(outmultinest):
		print("ERROR: The output dir: " + outmultinest + " does not exist!")
		sys.exit(1)
	if not os.path.isfile(outmultinest + config_name):
		shutil.copy(config_file, outmultinest + config_name)
	
	covfile = config['params']['covfile']

	if not os.path.isfile(covfile):
		print("ERROR: The cov file: " + covfile + " does not exist!")
		sys.exit(1)

	
	# LOAD the necessary data
	print('STATUS: Select the halos and read the reference 2PCF...')
	data_instance = DataClass(config_file)
	ref_data_vector = np.append(data_instance.ref_data_vector, data_instance.ndens_avg)
	print('STATUS: Loading the covariance matrix and invert it...')
	cov = np.loadtxt(covfile)
	icov = np.linalg.inv(cov)
	print(ref_data_vector.shape, icov.shape[0])

	if len(ref_data_vector) != icov.shape[0]:
		print("ERROR: The number of data bins is different than the size of the Cov matrix!")
		sys.exit(1)
	
	print('STATUS: Loading the low resolution data and create the model...')
	model_var = ModelClass(config_file, data_instance.halo_files)

	compute_2pcf_instance = ComputeCFClass(config_file)

	compute_pspec_instance = None#ComputePKClass(config_file)
	
	chi2_var = Chi2Class(config_file, model_var, ref_data_vector, icov, compute_2pcf_instance, compute_pspec_instance)

	#### BEGIN TEST
	#ic = 0.8 #np.zeros(1)
	#print("TEST_INFO: My test chi2 is: ", chi2_var.compute_chi2(ic))
	#sys.exit()
	#### END TEST

	minimizer_var = ClassMinimizer(config_file, chi2_var)
	#minimizer_var.one_param_minuit()
	#minimizer_var.five_params_minuit()
	#minimizer_var.nine_params_minuit()

	minimizer_var.scipy_dual_annealing()
	# minimizer_var.analyse_multinest()
	# minimizer_var.sample_values()
	#minimizer_var.scipy_minimizer()

	end = time.time()
	print('END: The code has finished in {time} seconds'.format(time=end - start))


if __name__ == '__main__':
	main()
