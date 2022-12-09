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
from apply_hod_cat import ModelClass


#python main_cat.py --case 1 1> out.case1.out 2> err.case1.err &
#python main_cat.py --case 2 1> out.case2.out 2> err.case2.err &
#python main_cat.py --case 5 1> out.case5.out 2> err.case5.err &
#python main_cat.py --case 6 1> out.case6.out 2> err.case6.err &

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


	# if case_ == 1:
		# config_file = f'{main_path}/config_RSD_diff_1296_1000_0.5_6pars_cf_7.5_47.5'
		# print(config_file)
		# best_fit_params = np.array([1.260433801744823690e+01, 2.080501379302197140e+00, 1.237742541378510452e+01, 3.393159184173180076e+00, 8.082312018031215795e-01, 9.328199250183341062e-01])

	# elif case_ == 2:
		# config_file = f'{main_path}/config_RSD_diff_1296_1000_0.5_6pars_cf_12.5_47.5'
		# print(config_file)
		# best_fit_params = np.array([1.219966361033310775e+01, 1.255253091182735758e+00, 1.179931479360136137e+01, 1.996987123804557385e+01, 6.600687949076093908e-01, 9.619679116420107867e-01])

	if case_ == 3:		
		config_file = f'{main_path}/config_RSD_diff_1536_1000_0.5_6pars_cf_7.5_47.5'
		print(config_file)
		best_fit_params = np.array([1.261826667540075775e+01, 2.007229157807212339e+00, 1.226249211011251106e+01, 2.676145552129969296e+00, 7.225904314671200979e-01, 9.687099501684696135e-01])

	elif case_ == 4:
		config_file = f'{main_path}/config_RSD_diff_1536_1000_0.5_6pars_cf_12.5_47.5'
		best_fit_params = np.array([1.244399906328524708e+01, 1.709861490971135378e+00, 1.250632813517288078e+01, 5.052549452135196795e+00, 9.005585122038467505e-01, 9.323787322622658635e-01])
		print(config_file)
	# elif case_ == 5:
		# config_file = f'{main_path}/config_RSD_diff_1296_1000_0.5_6pars_pk_0.02_0.3'
		# print(config_file)
		# best_fit_params = np.array([1.275575134641997721e+01, 2.302131016704958899e+00, 1.111421973396371676e+01, 2.194602840070038141e+00, 3.762552947398257186e-01, 1.042605643693187423e+00])

	# elif case_ == 6:
		# config_file = f'{main_path}/config_RSD_diff_1296_1000_0.5_6pars_pk_0.02_0.4'
		# print(config_file)
		# best_fit_params = np.array([1.263258316035693340e+01, 2.109388666996832917e+00, 1.219477256653407693e+01, 2.840436396269645503e+00, 6.860994928694214012e-01, 9.914476263000456813e-01])

	elif case_ == 7:
		config_file = f'{main_path}/config_RSD_diff_1536_1000_0.5_6pars_pk_0.02_0.3'
		best_fit_params = np.array([1.230054855283565196e+01, 1.401213148265698427e+00, 9.487868326362713489e+00, 1.234243376817737570e+01, 2.711188089205858187e-01, 1.048799591823725397e+00])
		print(config_file)
	elif case_ == 8:
		config_file = f'{main_path}/config_RSD_diff_1536_1000_0.5_6pars_pk_0.02_0.4'
		print(config_file)
		best_fit_params = np.array([1.245062848525131649e+01, 1.690872606626584895e+00, 1.101067240543193648e+01, 6.041102985260993208e+00, 3.962846331574765291e-01, 1.031862829655151526e+00])


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
	

	halo_file_pattern = config['params']['halofile']
	halo_files = glob.glob(halo_file_pattern.format("*"))
	print(halo_file_pattern)
	print(len(halo_files))
	outpath = config['params']["outpathcat"]

	

	
	# LOS = [902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921]
	# halo_files = []
	# for LOS_ in LOS:
		# halo_files.append(halo_file_pattern.format(LOS_))

	halo_files = glob.glob(halo_file_pattern.format("*"))	
	print(len(halo_files))

	print('STATUS: Loading the low resolution data and create the model...')


	print(len(halo_files))
	for i in range(len(halo_files)):
		filename = os.path.basename(halo_files[i: i + 1][0])

		if not os.path.isfile(outpath  + filename):
			model_var = ModelClass(config_file, halo_files[i: i + 1])
			dict_of_gsamples, dict_gal_type = model_var.populate_mock(best_fit_params, 0.00421339953605372, i)
			

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

