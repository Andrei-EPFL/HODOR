import sys
import numpy as np
# from halotools.mock_observables import tpcf_multipole
# from halotools.mock_observables import s_mu_tpcf
# import Corrfunc.theory as cft


sys.path.append("/global/homes/a/avariu/phd/danielscodes/pyfcfc/")
from pyfcfc.boxes import py_compute_cf

class ComputeCFClass():
    """ The class that computes the 2PCF given the config file and a sample of galaxies """
    def __init__(self):

        self.s_min = 0
        self.s_max = 50
        self.n_s_bins = 10

        self.s_bins = np.linspace(self.s_min, self.s_max, self.n_s_bins + 1).astype(np.double)
        print("INFO: The s_bins are: ", self.s_bins)
        self.r = np.array([(a + b) / 2.0 for a, b in zip(self.s_bins[:-1], self.s_bins[1:])])

        self.mu_min = 0
        self.mu_max = 1
        self.n_mu_bins = 100
        self.mu_bins = np.linspace(self.mu_min, self.mu_max, self.n_mu_bins + 1)
        print("INFO: The mu_bins are: ", self.mu_bins)

        self.nthreads = 64
        self.box_size = 505


    def compute_2pcf_pyfcfc(self, x_c, y_c, z_c):
        """ Compute the 2PCF of galaxies given the
        dictionary of galaxies. It uses the corrfunc package (which is much faster than halotools)
        to compute the 2PCF, with a natural estimator """
        
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
        
        xi0 = np.array(results["multipoles"][0][0]).copy()
        xi2 = np.array(results["multipoles"][0][1]).copy()
        xi4 = np.array(results["multipoles"][0][2]).copy()

        return self.r, xi0, xi2, xi4


    def compute_2pcf_halotools(self, x_c, y_c, z_c):
        """ Compute the 2PCF of galaxies given the
        dictionary of galaxies. It uses the halotools package to compute the 2PCF """

        sample1 = np.array([x_c, y_c, z_c]).T
        print(sample1.shape)

        print(self.s_bins)
        # Compute xi0, xi2, xi4
        xi_s_mu = s_mu_tpcf(sample1, self.s_bins, self.mu_bins, period=self.box_size, num_threads=self.nthreads)
        xi0 = tpcf_multipole(xi_s_mu, self.mu_bins, order=0)
        xi2 = tpcf_multipole(xi_s_mu, self.mu_bins, order=2)
        xi4 = tpcf_multipole(xi_s_mu, self.mu_bins, order=4)

        return self.r, xi0, xi2, xi4


    def spherical_sector_volume(self, s, mu):
        """
        This function is used to calculate analytical randoms.
        Calculate the volume of a spherical sector, used for the analytical randoms.
        https://en.wikipedia.org/wiki/Spherical_sector
        Note that the extra factor of 2 is to get the reflection.
        """
        theta = np.arccos(mu)

        vol = (2.0*np.pi/3.0) * np.outer((s**3.0), (1.0-np.cos(theta)))*2.0
        return vol


    def compute_analytical_RR(self, NR):
        """ This function is used to compute an analytical value of the
        RR term of the 2PCF """

        # do volume calculations
        mu_bins_reverse_sorted = np.sort(self.mu_bins)[::-1]
        dv = self.spherical_sector_volume(self.s_bins, mu_bins_reverse_sorted)
        dv = np.diff(dv, axis=1)  # volume of wedges
        dv = np.diff(dv, axis=0)  # volume of wedge 'pieces'
        global_volume = self.box_size**3

        # calculate the random-random pairs.
        rhor = NR**2/global_volume
        RR = (dv*rhor)
        return RR


    def compute_2pcf_corrfunc(self, x_c, y_c, z_c):
        """ Compute the 2PCF of galaxies given the
        dictionary of galaxies. It uses the corrfunc package (which is much faster than halotools)
        to compute the 2PCF, with a natural estimator """

        RR = self.compute_analytical_RR(len(x_c))

        autocorr = 1
        DD_s_mu = cft.DDsmu(autocorr, self.nthreads, self.s_bins, self.mu_max, self.n_mu_bins, x_c, y_c, z_c, boxsize=self.box_size, verbose=True)
        DD_s_mu = np.reshape(DD_s_mu[:]['npairs'], (self.n_s_bins, self.n_mu_bins))
        # Corrfunc adds also the central galaxy in an s bin (when smin==0). So one needs to subtract it from the number of pairs.
        if self.s_min == 0:
            DD_s_mu[0][0] = DD_s_mu[0][0] - len(x_c)
            
        xi_s_mu = (DD_s_mu / RR) - 1

        xi0 = tpcf_multipole(xi_s_mu, self.mu_bins, order=0)
        xi2 = tpcf_multipole(xi_s_mu, self.mu_bins, order=2)
        xi4 = tpcf_multipole(xi_s_mu, self.mu_bins, order=4)

        return self.r, xi0, xi2, xi4
    

def main():
    import glob
    import os

    cf_inst =  ComputeCFClass()
    path = "/global/homes/a/avariu/desi_project_dir/FastPM_SLICS/SLICS/galaxy/*"
    
    files = glob.glob(path)
    for file_ in files:
        x_c, y_c, z_c = np.loadtxt(file_, usecols=(0, 1, 2), unpack=True)
        s, cf0, cf2, cf4 = cf_inst.compute_2pcf_pyfcfc(x_c, y_c, z_c)
        np.savetxt(f"/global/homes/a/avariu/desi_project_dir/FastPM_SLICS/SLICS/2pcf/real/{os.path.basename(file_)}", np.array([s, cf0, cf2, cf4]).T)

    # s, cf0, cf2, cf4 = cf_inst.compute_2pcf_corrfunc(x_c, y_c, z_c)
    # np.savetxt("./test/1.041halo.dat_LOS985_corrfunc.gcat", np.array([s, cf0, cf2, cf4]).T)

def plot_aux(ax, file, usecols, color, label):
    s, cf0, cf2, cf4 = np.loadtxt(file, usecols=usecols, unpack=True)

    range_ = (s >=0.0) & (s <=50)
    s, cf0, cf2, cf4  = s[range_], cf0[range_], cf2[range_], cf4[range_]

    ax[0].plot(s, s*s* cf0, color=color, label=label)
    ax[1].plot(s, s*s* cf2, color=color, label=label)
    ax[2].plot(s, s*s* cf4, color=color, label=label)

    return s, cf0, cf2, cf4

def plot_aux_err(ax, file, color, label):
    s, cf0, cf2, cf4, std0, std2, std4 = np.loadtxt(file, usecols=(0, 1, 2, 3, 4, 5, 6), unpack=True)
    range_ = (s >=0.0) & (s <=50)
    s, cf0, cf2, cf4, std0, std2, std4  = s[range_], cf0[range_], cf2[range_], cf4[range_], std0[range_], std2[range_], std4[range_]

    ax[0].errorbar(s, s*s*cf0, yerr=s*s*std0, color=color, label=label)
    ax[1].errorbar(s, s*s*cf2, yerr=s*s*std2, color=color, label=label)
    ax[2].errorbar(s, s*s*cf4, yerr=s*s*std4, color=color, label=label)

    return std0* np.sqrt((0.5 ** 3) / 24) , std2* np.sqrt((0.5 ** 3) / 24) , std4* np.sqrt((0.5 ** 3) / 24) 

def plot():
    import matplotlib.pyplot as pt

    fig, ax = pt.subplots(2, 3, figsize=(15, 10), sharex=True)

    so, cfr0, cfr2, cfr4 = plot_aux(ax[0], "/global/project/projectdirs/desi/cosmosim/cov_mat_challenge/cubic_box/SLICS/CF_multipoles_V2/rsd/xil_1.041halo_LOS985.txt", (0,1,2,3), "orange", "SLICS_V2")

    std0, std2, std4 = plot_aux_err(ax[0], "/global/homes/a/avariu/desi_project_dir/FastPM_SLICS/scripts/avg_std_slics.dat", color="black", label="avg_slics_v2")

    sr, cfo0, cfo2, cfo4 = plot_aux(ax[0], "./test/1.041halo.dat_LOS985_corrfunc.gcat", (0, 1, 2, 3), "red", "corrfunc")
    ax[1][0].plot(sr, (cfr0 - cfo0)/std0, color="red")
    ax[1][1].plot(sr, (cfr2 - cfo2)/std2, color="red")
    ax[1][2].plot(sr, (cfr4 - cfo4)/std4, color="red")

    sn, cfn0, cfn2, cfn4 = plot_aux(ax[0], "./test/1.041halo.dat_LOS985_pyfcfc.gcat", (0, 1, 2, 3), "blue", "pyfcfc")
    ax[1][0].plot(sr, (cfr0 - cfn0)/std0, color="blue", ls="--")
    ax[1][1].plot(sr, (cfr2 - cfn2)/std2, color="blue", ls="--")
    ax[1][2].plot(sr, (cfr4 - cfn4)/std4, color="blue", ls="--")

    sc, cfc0, cfc2, cfc4 = plot_aux(ax[0], "/global/project/projectdirs/desi/cosmosim/cov_mat_challenge/cubic_box/SLICS/CF_multipoles_V2_clean/rsd/xil_1.041halo_LOS985.txt", (0,1,2,3), "green", "SLICS_V2_clean")
    ax[1][0].plot(sr, (cfr0 - cfc0)/std0, color="orange")
    ax[1][1].plot(sr, (cfr2 - cfc2)/std2, color="orange")
    ax[1][2].plot(sr, (cfr4 - cfc4)/std4, color="orange")


    # ax[1][1].set_ylim([-0.1, 0.1])
    ax[1][0].set_ylabel("s*s*(difference)")

    ax[0][0].legend()
    fig.savefig("./test/tpcf.png")

if __name__ == '__main__':
    main()
    # plot()
