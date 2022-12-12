""" Apply HOD and fit the 2PCF to High resolution simulation """
#!/global/homes/a/avariu/.conda/envs/hodfit/bin/python
import sys
import os
import shutil
import glob
import time
import argparse
import configparser
import numpy as np
from hod_pack import ModelClass


def main():
	""" Main function of the code"""
	# Parameters
	start = time.time()
	parser = argparse.ArgumentParser()
	parser.add_argument("--config", type=str, help="ini file holding configuration")
	parser.add_argument("--case", type=int, help="1 to 8")

	args = parser.parse_args()
	case_ = args.case # config file

	#### Configuration file
	main_path = "/global/homes/a/avariu/phd/hodsimulations/HOD_pipeline_Cheng/fit/configs/FastPM4SLICS/for_cat/"

	config_file = f'{main_path}/config_RSD_diff_1536_1000_0.5_6pars_pk_0.02_0.4'
	print(config_file)
	best_fit_params = np.array([1.245062848525131649e+01, 1.690872606626584895e+00, 1.101067240543193648e+01, 6.041102985260993208e+00, 3.962846331574765291e-01, 1.031862829655151526e+00])

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
	

	halo_file_pattern = config['params']['halo_file']
	halo_files = glob.glob(halo_file_pattern.format("*"))
	print(halo_file_pattern)
	print(len(halo_files))
	outpath = config['params']["outpath_catalog"]


	# LOS = [902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921]
	# halo_files = []
	# for LOS_ in LOS:
		# halo_files.append(halo_file_pattern.format(LOS_))

	# halo_files = glob.glob(halo_file_pattern.format("*"))	
	# print(len(halo_files))

	print('STATUS: Loading the low resolution data and create the model...')


	print(len(halo_files))
	for i in range(len(halo_files)):
		filename = os.path.basename(halo_files[i: i + 1][0])

		if not os.path.isfile(outpath  + filename):
			model_var = ModelClass(config_file, halo_files[i: i + 1])
			dict_of_gsamples = model_var.populate_mock(best_fit_params, 0.00421339953605372, indx=i)
			

			for j, key in enumerate(dict_of_gsamples.keys()):
				print(halo_files[i: i + 1], i, key, flush=True)
				np.savetxt(outpath  + filename, dict_of_gsamples[key])
			
	param_names = ' '.join(model_var.parameters_names)
	np.savetxt(outpath + "/best_fit_params.txt", best_fit_params.reshape(1, best_fit_params.size), header=param_names)
	
	# 	# for j, key in enumerate(dict_of_gsamples.keys()):

	# 		# outfile = outpath + key[0: -10]  + "_gal.txt"
	# 		# print(j, key, outfile)
	# 	pylians_pk = compute_pspec_instance.pspec_pylians(dict_of_gsamples, i)
	# 	# nbody_pk = compute_pspec_instance.pspec_2d(dict_of_gsamples, i)
		
	# 	# outfile_n = outpath + "/nbody_pk.txt"
	# 	outfile_p = outpath + "/pylians_pk_cic.txt"
		

	# 	# np.savetxt(outfile_n, nbody_pk) 
	# 	np.savetxt(outfile_p, pylians_pk) 
	
	# 	exit()
	
	end = time.time()
	print('END: The code has finished in {time} seconds'.format(time=end - start))


if __name__ == '__main__':
	main()

