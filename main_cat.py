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
from hod_pack import ModelClass, DataClass


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

	config_file = "./configs/config.ini"
	print(config_file)
	# best_fit_params = np.array([1.247878283003146649e+01, 1.780774210878682151e+00, 1.114040094979609385e+01, 4.845500406517473380e+00, 3.915263394920053264e-01, 1.094422425428754231e+00])
	best_fit_params = np.array([12.73944303708267, 2.2775967871074885, 11.999032351426273, 1.880455318916949, 0.5778940798839721, 1.017176352435166])
	
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
	# halo_files = glob.glob(halo_file_pattern.format("*"))
	print(halo_file_pattern)
	# print(len(halo_files))
	outpath = config['params']["outpath_catalog"]

	LOS_list = [902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 923, 925, 926, 954, 955, 957, 959, 960, 961, 964, 965, 966, 967, 968, 969, 970, 971, 972, 973, 974, 977, 978, 979, 980, 981, 982, 983, 984, 985, 986, 988, 989, 990, 992, 993, 994, 996, 997, 998, 999, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1014, 1015, 1017, 1019, 1022, 1024, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1038, 1044, 1045, 1047, 1048, 1049, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1066, 1070, 1071, 1072, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1098, 1099, 1100]

	# LOS_list = [902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921]
	halo_files = []
	for LOS_ in LOS_list:
		halo_files.append(halo_file_pattern.format(LOS_))

	# halo_files = glob.glob(halo_file_pattern.format("*"))	
	# print(len(halo_files))

	print('STATUS: Loading the low resolution data and create the model...')

	# start = time.time()
	# data_instance = DataClass(config_file, halo_files=["/global/cscratch1/sd/avariu/fastpm_nonfixed-amplitude_L2000_N5184_S537-B2-T40_fof_0.2500.hdf5"])
	# print(time.time() - start, flush=True)

	# model_var = ModelClass(config_file, data_instance.halo_files, data_instance.halo_cat_list)
	# print(time.time() - start, flush=True)

	# dict_of_gsamples = model_var.populate_mock(best_fit_params, 0.00421339953605372, indx=1)
	# for j, key in enumerate(dict_of_gsamples.keys()):
	# 	print(dict_of_gsamples[key].shape)
	# print(time.time() - start, flush=True)
	# exit()

	# print(len(halo_files))
	# start = time.time()
	
	for i in range(len(halo_files)):
		filename = os.path.basename(halo_files[i: i + 1][0])
		print(filename)
		if not os.path.isfile(outpath  + filename):
			data_instance = DataClass(config_file, halo_files=halo_files[i: i + 1])

			model_var = ModelClass(config_file, data_instance.halo_files, data_instance.halo_cat_list)
	
			dict_of_gsamples = model_var.populate_mock(best_fit_params, 0.00421339953605372, indx=i)

			for j, key in enumerate(dict_of_gsamples.keys()):
					print(dict_of_gsamples[key].shape)
					np.savetxt(outpath  + filename, dict_of_gsamples[key], fmt='%f %f %f %f %d')
	
	# print(time.time() - start, flush=True)

	# print(len(halo_files))
	# start = time.time()
	# data_instance = DataClass(config_file, halo_files=halo_files)
	# print(time.time() - start, flush=True)

	# model_var = ModelClass(config_file, data_instance.halo_files, data_instance.halo_cat_list)
	# dict_of_gsamples = model_var.populate_mock(best_fit_params, 0.00421339953605372, indx=1)
			
	# for j, key in enumerate(dict_of_gsamples.keys()):
	# 	print(dict_of_gsamples[key].shape)
	# 	np.savetxt(outpath  + filename, dict_of_gsamples[key], fmt='%f %f %f %f %d')
			
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


