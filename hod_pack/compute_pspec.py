import sys
import configparser
import numpy as np

sys.path.append("/global/homes/a/avariu/phd/danielscodes/pypowspec/")
from pypowspec import compute_auto_box

class ComputePKClass():
    """ The class that computes the pspec given the config file and a sample of galaxies """
    def __init__(self, config_file):
        config = configparser.ConfigParser()
        config.read(config_file)

        self.fit_kmin = config['pspec'].getfloat('fit_kmin')
        self.fit_kmax = config['pspec'].getfloat('fit_kmax')
        self.box_size = config['params'].getfloat('box_size')

        self.config_file = config['pspec']["config_file"]
        if not os.path.isfile(self.config_file):
            print("ERROR: The power spectrum config file is not correctly give", flush=True)
            sys.exit(1)

    def pspec(self, dict_of_gsamples, index):
        ### Multipoles list
        pk0_list = []
        pk2_list = []
        pk4_list = []
        
        ### Go through all galaxy samples
        for j, key in enumerate(dict_of_gsamples.keys()):
            x_c, y_c, z_c = dict_of_gsamples[key]["x"], dict_of_gsamples[key]["y"], dict_of_gsamples[key]["z_rsd"]
            
            x_c = (x_c + self.box_size) % self.box_size
            y_c = (y_c + self.box_size) % self.box_size
            z_c = (z_c + self.box_size) % self.box_size
            
            w = np.ones(len(x_c))
            pk = compute_auto_box(x_c.astype(np.float), y_c.astype(np.float), z_c.astype(np.float), w.astype(np.float), powspec_conf_file=self.config_file)
        
            k = pk["k"]
            pk0 = pk["multipoles"][:, 0]
            pk2 = pk["multipoles"][:, 1]
            pk4 = pk["multipoles"][:, 2]

            pk0_list.append(pk0)
            pk2_list.append(pk2)
            pk4_list.append(pk4)

        range_ = (self.fit_kmin <= k) & (k <= self.fit_kmax)

        mean_pk0 = np.mean(np.array(pk0_list), 0)[range_]
        mean_pk2 = np.mean(np.array(pk2_list), 0)[range_]
        mean_pk4 = np.mean(np.array(pk4_list), 0)[range_]

        return k[range_], mean_pk0, mean_pk2, mean_pk4

