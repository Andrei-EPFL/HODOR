import configparser
import numpy as np
# import scipy.optimize as opt
# from iminuit import Minuit
import pymultinest as pmn

class ClassMinimizer():
	""" A class containing more minimizers """
	def __init__(self, config_file, chi2_instance):
		config = configparser.ConfigParser()
		config.read(config_file)

		self.chi2_instance = chi2_instance
		self.ic = np.array([+10, 4, 12, 3, 0.5])#np.zeros(self.num_params) #[+0.4, 0, -1, -0.1, 0.71]#

		self.num_sampled_values = config['params'].getint('num_sampled_values')

		self.num_params        = config['params'].getint('num_params')
		self.outpath_info      = config['params']['outpath_info']
		self.outpath_multinest = config['params']['outpath_multinest']

		self.live_points    = config["multinest"].getint("live_points")
		self.tol            = config["multinest"].getfloat("tol")
		self.verbose        = config["multinest"].getboolean("verbose")

		self.parameters_names = self.chi2_instance.parameters_names
		self.prior_params = {}
		for p in self.parameters_names:
			self.prior_params[p] = tuple(map(float, config.get('priors', p).split(', ')))
		print("INFO: The priors of the parameters are " + str(self.prior_params))

		partname = ""
		for index in range(self.num_params):
			l = self.prior_params[self.parameters_names[index]][0]
			u = self.prior_params[self.parameters_names[index]][1]
			partname = partname + self.parameters_names[index] + "_" + str(l) + "_" + str(u) + "_"

		self.outbase = self.outpath_multinest + "MN_" + partname
		print("INFO: The output for multinest is ", self.outbase)


	def write_info_minuit(self, m):
		""" Write the info provided by Minuit """
		file = open(self.outpath_info+'/Minuit_results.txt', 'w')
		file.write(str(m.get_fmin())+'\n')
		file.write(str(m.values)+'\n')
		file.write(str(m.errors)+'\n')
		file.write('0.5*chi2 : '+str(m.fval)+'\n')
		file.write(str(m.get_param_states())+'\n\n\n')

		content = '#chi2_bestfit: ' + str(2*m.fval) + "\n"
		for i in range(self.num_params):
			content = content+' '+ m.get_param_states()[i].name +"=" + str(m.np_values()[i])+' errh_'+m.get_param_states()[i].name+ "=" + str(m.np_errors()[i])
			content += ' \n'

		file.write(content)
		file.close()


	def one_param_minuit(self):
		""" Minuit: For one parameter, use this """
		l = self.prior_params[self.parameters_names[0]][0]
		u = self.prior_params[self.parameters_names[0]][1]
		print("STATUS: Start minimizing...")
		m = Minuit(self.chi2_instance.compute_chi2, forced_parameters=('a'), a=self.ic[0], \
			error_a=0.1, limit_a=((l, u)), print_level=1, errordef=0.5)

		m.migrad()

		self.write_info_minuit(m)


	def nine_params_minuit(self):
		""" Minuit: For 9 parameters, use this """
		l = np.ones(self.num_params)
		u = np.ones(self.num_params)
		for i in range(self.num_params):
			l[i] = self.prior_params[self.parameters_names[i]][0]
			u[i] = self.prior_params[self.parameters_names[i]][1]

		print("STATUS: Start minimizing...")
		m = Minuit.from_array_func(self.chi2_instance.compute_chi2, self.ic, \
			error=(np.zeros(self.num_params) + 0.1), \
			fix=(False, False, False, False, False, False, False, False, False), \
			limit=((l[0], u[0]), (l[1], u[1]), (l[2], u[2]), (l[3], u[3]), (l[4], u[4]), (l[5], u[5]), (l[6], u[6]), (l[7], u[7]), (l[8], u[8])), \
			name=("Fa", "Fb", "logMc", "siga", "sigb", "Fs", "logMmin", "deltaM", "alpha"), \
			print_level=1, errordef=1)

		m.migrad()

		self.write_info_minuit(m)


	def five_params_minuit(self):
		""" Minuit: For 5 parameters, use this """
		l = np.ones(self.num_params)
		u = np.ones(self.num_params)
		for i in range(self.num_params):
			l[i] = self.prior_params[self.parameters_names[i]][0]
			u[i] = self.prior_params[self.parameters_names[i]][1]

		print("STATUS: Start minimizing...")
		m = Minuit.from_array_func(self.chi2_instance.compute_chi2, self.ic, \
			error=(np.zeros(self.num_params) + 0.1), \
			fix=(False, False, False, False, False, False, False, False, False), \
			limit=((l[0], u[0]), (l[1], u[1]), (l[2], u[2]), (l[3], u[3]), (l[4], u[4])), \
			name=("logMcut", "sigma_logM", "logM1", "k", "alpha"), \
			print_level=1, errordef=1)

		m.migrad()

		self.write_info_minuit(m)


	def scipy_minimizer(self):
		""" Scipy minimizer"""
		l = np.ones(self.num_params)
		u = np.ones(self.num_params)
		for i in range(self.num_params):
			l[i] = self.prior_params[self.parameters_names[i]][0]
			u[i] = self.prior_params[self.parameters_names[i]][1]
		bounds_ = opt.Bounds(l, u)

		print("STATUS: Start minimizing...")
		res = opt.minimize(self.chi2_instance.compute_chi2, self.ic, tol=1e-3, bounds=bounds_, method='Powell')
		best_fit_params = res.x
		print('INFO: Best-fit param: ', best_fit_params)
		print('with a minimum chi2 %f ' % (self.chi2_instance.compute_chi2(best_fit_params)))
		print('INFO: Success:', res.success)
		print('INFO: Message:', res.message)


	def scipy_dual_annealing(self):
		""" Scipy dual annealing """
		l = np.ones(self.num_params)
		u = np.ones(self.num_params)
		for i in range(self.num_params):
			l[i] = self.prior_params[self.parameters_names[i]][0]
			u[i] = self.prior_params[self.parameters_names[i]][1]
		# bounds_ = opt.Bounds(l, u)
		print(np.vstack((l, u)).T.shape)

		res = opt.dual_annealing(self.chi2_instance.compute_chi2, np.vstack((l, u)).T)
		best_fit_params = res.x
		print('INFO: Best-fit param: ', best_fit_params)
		print('with a minimum chi2 %f ' % (self.chi2_instance.compute_chi2(best_fit_params)))
		print('INFO: Success:', res.success)
		print('INFO: Message:', res.message)


	def sample_values(self):
		""" simple sampling of chi2 with values between min and max
		for ONE parameter """
		index = 3
		l = self.prior_params[self.parameters_names[index]][0]
		u = self.prior_params[self.parameters_names[index]][1]
		# partname = "best_fitrsdspL0L2smin12.5_multinest_30_rsd_1536_vdisp"
		values = np.linspace(l, u, num=self.num_sampled_values)
		print('INFO: The sampled values are : ', values)
		chi2_arr = np.zeros(values.shape)
		nsat_arr = np.zeros(values.shape)
		ncen_arr = np.zeros(values.shape)

		for i, val_ in enumerate(values):
			
			partname = "1296_1_test"

			chi2_arr[i], nsat_arr[i], ncen_arr[i] = self.chi2_instance.return_chi2_save_clustering([1.244424312095842922e+01, 3.780612391087811108e+00, 1.259489112661007049e+01, 3.662500985959312061e+00, 4.096474641170874120e-01], partname)
			# self.chi2_instance.save_galaxy_catalogs([1.244424312095842922e+01, 3.780612391087811108e+00, 1.259489112661007049e+01, 3.662500985959312061e+00, 9.096474641170874120e-01], partname)
	
			# chi2_arr[i], nsat_arr[i], ncen_arr[i] = self.chi2_instance.return_chi2_save_clustering([0.142785750296551424E+02,    0.484354250890535365E+01,    0.128656562324185906E+02,    0.146517358593474949E+02,    0.647165718765783771E+00], partname)
			# self.chi2_instance.save_galaxy_catalogs([0.142785750296551424E+02,    0.484354250890535365E+01,    0.128656562324185906E+02,    0.146517358593474949E+02,    0.647165718765783771E+00], partname)
			
		ratio_arr = nsat_arr / (nsat_arr + ncen_arr)
		np.savetxt(self.outpath_info + "/Sample_val_chi2_" + partname + "_.txt", np.array([1, chi2_arr, nsat_arr, ncen_arr, ratio_arr]).T)


	def prior(self, cube, ndim, nparams):
		for i in range(self.num_params):
			val_max = self.prior_params[self.parameters_names[i]][1]
			val_min = self.prior_params[self.parameters_names[i]][0]
			cube[i] = cube[i] * (val_max - val_min) + val_min


	def loglike(self, cube, ndim, nparams):
		lnlike = -0.5 * self.chi2_instance.compute_chi2(cube[0:self.num_params])
		return lnlike


	def run_multinest(self):
		print("\nSTATUS: Begin Multinest...")
		pmn.run(self.loglike, self.prior, self.num_params, outputfiles_basename=self.outbase, resume=True, \
			verbose=self.verbose, n_live_points=self.live_points, evidence_tolerance=self.tol, seed=79)


	def analyse_multinest(self):
		res = pmn.Analyzer(outputfiles_basename=self.outbase, n_params=self.num_params)
		best_fit_params = res.get_best_fit()['parameters']

		[s, xi0, xi2, xi4], [k, pk0, pk2, pk4], chi2, [nsat, ncen] = self.chi2_instance.compute_best_fit(best_fit_params)

		print(f"INFO: The chi2 for {self.outbase} is equal to {chi2} ")
		print(("INFO: The best fit parameters are: [" + ', '.join(['%f']*len(best_fit_params))+"]") % tuple(best_fit_params))

		best_pars_str = ' '.join(['%f']*len(best_fit_params)) % tuple(best_fit_params)  + f"; chi2={chi2}"

		np.savetxt(self.outbase + 'best_cf.dat', np.transpose([s, xi0, xi2, xi4]), fmt='%.8g', header=best_pars_str)
		np.savetxt(self.outbase + 'best_pk.dat', np.transpose([k, pk0, pk2, pk4]), fmt='%.8g', header=best_pars_str)

		ratio_arr = nsat / (nsat + ncen)
		np.savetxt(self.outpath_info + "/chi2_nsat_ncen_satfrac.txt", np.array([chi2, nsat, ncen, ratio]).T)

		param_names = ' '.join(self.parameters_names)
		np.savetxt(self.outpath_info + "/best_fit_params.txt", np.array(best_fit_params).reshape(1, len(best_fit_params)), header=param_names)
