import sys
import configparser
import numpy as np

sys.path.append("/global/homes/a/avariu/phd/danielscodes/pyfcfc/")
from pyfcfc.boxes import py_compute_cf

class ComputeCFClass():
    """ The class that computes the 2PCF given the config file and a sample of galaxies """
    def __init__(self, config_file):
        config = configparser.ConfigParser()
        config.read(config_file)

        self.verbose = config['params'].getboolean('verbose')
        self.s_min = config['TPCF'].getfloat("s_min")
        self.s_max = config['TPCF'].getfloat("s_max")
        self.n_s_bins = config['TPCF'].getint("n_s_bins")
        #self.s_bins = np.geomspace(self.s_min, self.s_max, self.n_s_bins + 1)
        self.s_bins = np.linspace(self.s_min, self.s_max, self.n_s_bins + 1).astype(np.double)
        print("INFO: The s_bins are: ", self.s_bins)
        self.r = np.array([(a + b) / 2.0 for a, b in zip(self.s_bins[:-1], self.s_bins[1:])])

        self.mu_min = config['TPCF'].getfloat("mu_min")
        self.mu_max = config['TPCF'].getfloat("mu_max")
        self.n_mu_bins = config['TPCF'].getint("n_mu_bins")
        self.mu_bins = np.linspace(self.mu_min, self.mu_max, self.n_mu_bins + 1)
        print("INFO: The mu_bins are: ", self.mu_bins)

        self.nthreads = config['TPCF'].getint("nthreads")
        self.box_size = config['params'].getfloat('box_size')


    def compute_2pcf_pyfcfc(self, dict_of_gsamples, index):
        """ Compute the 2PCF of galaxies given the
        dictionary of galaxies."""

        xi0_arr = np.zeros((len(dict_of_gsamples.keys()), len(self.r)))
        xi2_arr = np.zeros((len(dict_of_gsamples.keys()), len(self.r)))
        xi4_arr = np.zeros((len(dict_of_gsamples.keys()), len(self.r)))

        for j, key in enumerate(dict_of_gsamples.keys()):            
            x_c, y_c, z_c = dict_of_gsamples[key]["x"], dict_of_gsamples[key]["y"], dict_of_gsamples[key]["z_rsd"]
            x_c = (x_c + self.box_size) % self.box_size
            y_c = (y_c + self.box_size) % self.box_size
            z_c = (z_c + self.box_size) % self.box_size

            xyz_gal = np.array([x_c, y_c, z_c]).T.astype(np.double)
            print(xyz_gal.shape)
            # Compute xi0, xi2, xi4
            results = py_compute_cf([xyz_gal], [np.ones(xyz_gal.shape[0])], 
                        self.s_bins.copy(), 
                        None, 
                        self.n_mu_bins, 
                        label = ['D'], # Catalog labels matching the number of catalogs provided
                        bin=1, # bin type for multipoles
                        pair = ['DD'], # Desired pair counts
                        box=self.box_size,
                        multipole = [0, 2, 4], # Multipoles to compute
                        cf = ['DD / @@ - 1'])			
            
            print(results["s"])
            print(np.array(results["multipoles"][0][0][:]))
            xi0_arr[j] = np.array(results["multipoles"][0][0]).copy()
            xi2_arr[j] = np.array(results["multipoles"][0][1]).copy()
            xi4_arr[j] = np.array(results["multipoles"][0][2]).copy()

        mean_xi0 = np.mean(xi0_arr, 0)
        mean_xi2 = np.mean(xi2_arr, 0)
        mean_xi4 = np.mean(xi4_arr, 0)
        return self.r, mean_xi0, mean_xi2, mean_xi4
    