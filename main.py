""" Apply HOD and fit the 2PCF to High resolution simulation """
#!/global/homes/a/avariu/.conda/envs/hodfit/bin/python
import sys
import os
import shutil
import time
import argparse
import configparser
import numpy as np

from hod_pack import *

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
		config_file = './configs/config.ini'

	if not os.path.isfile(config_file):
		print(f"ERROR: The config file: {config_file} does not exist!")
		sys.exit(1)

	config_name = os.path.basename(config_file)

	config = configparser.ConfigParser()
	config.read(config_file)

	outpath_multinest = config['params']['outpath_multinest']
	if not os.path.isdir(outpath_multinest):
		print(f"ERROR: The output dir: {outpath_multinest} does not exist!")
		sys.exit(1)

	if not os.path.isfile(outpath_multinest + config_name):
		shutil.copy(config_file, outpath_multinest + config_name)

	covariance_file = config['params']['covariance_file']
	if not os.path.isfile(covariance_file):
		print(f"ERROR: The cov file: {covariance_file} does not exist!")
		sys.exit(1)

	# LOAD the necessary data
	print('STATUS: Select the halos and read the reference...')
	data_instance = DataClass(config_file)
	print('STATUS: Loading the covariance matrix and invert it...')
	cov = np.loadtxt(covariance_file)
	icov = np.linalg.inv(cov)
	print(data_instance.ref_data_vector.shape, icov.shape[0])

	if len(data_instance.ref_data_vector) != icov.shape[0]:
		print("ERROR: The number of data bins is different than the size of the Cov matrix!")
		sys.exit(1)
	
	dens_instance = DensClass(config_file, data_instance.halo_mass_files)

	print('STATUS: Create the model...')
	model_instance = ModelClass(config_file, data_instance.halo_files, data_instance.halo_cat_list)

	compute_2pcf_instance = ComputeCFClass(config_file)

	compute_pspec_instance = ComputePKClass(config_file)
	
	chi2_var = Chi2Class(config_file, model_instance, data_instance.ref_data_vector, icov, compute_2pcf_instance, compute_pspec_instance, dens_instance)

	#### BEGIN TEST
	### This test is conceived for the 5 parameter HOD model, together with the velocity dispersion
	# parameters = [1.247878283003146649e+01, 1.780774210878682151e+00, 1.114040094979609385e+01, 4.845500406517473380e+00, 3.915263394920053264e-01, 1.094422425428754231e+00]
	# print("TEST_INFO: My test chi2 is: ", chi2_var.compute_chi2(parameters))
	# sys.exit()
	#### END TEST

	minimizer_var = ClassMinimizer(config_file, chi2_var)
	#minimizer_var.one_param_minuit()
	#minimizer_var.five_params_minuit()
	#minimizer_var.nine_params_minuit()

	minimizer_var.run_multinest()
	# minimizer_var.analyse_multinest()
	# minimizer_var.sample_values()
	#minimizer_var.scipy_minimizer()

	end = time.time()
	print('END: The code has finished in {time} seconds'.format(time=end - start), flush=True)


if __name__ == '__main__':
	main()
