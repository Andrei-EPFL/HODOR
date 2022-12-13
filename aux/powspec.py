import sys
import numpy as np
from ctypes import *


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
    def __init__(self):

        self.box_size = 505


    def compute_catalog(self, x, y, z):

        x = x % self.box_size
        y = y % self.box_size
        z = z % self.box_size

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


    def pspec_old(self, x_c, y_c, z_c):
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

        
        ### Go through all galaxy samples
        cata = self.compute_catalog(x_c, y_c, z_c)

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
        return k, pk0, pk2, pk4


    def pspec_new(self, x_c, y_c, z_c, outfile):
        sys.path.append("/global/homes/a/avariu/phd/danielscodes/pypowspec/")
        from pypowspec import compute_auto_box

        x_c = (x_c + self.box_size) % self.box_size
        y_c = (y_c + self.box_size) % self.box_size
        z_c = (z_c + self.box_size) % self.box_size
        
        w = np.ones(len(x_c))

        pk = compute_auto_box(x_c, y_c, z_c, w, powspec_conf_file="/global/homes/a/avariu/desi_project_dir/FastPM_SLICS/config_files/pypowspec_hod.conf", output_file=f"./test/{outfile}")       

        k = pk["k"]
        pk0 = pk["multipoles"][:, 0]
        pk2 = pk["multipoles"][:, 1]
        pk4 = pk["multipoles"][:, 2]



def main():
    pk_inst =  ComputePKClass()
    path = "/global/homes/a/avariu/desi_project_dir/FastPM_SLICS/SLICS/galaxy/1.041halo.dat_LOS985.gcat"
    
    x_c, y_c, z_c = np.loadtxt(path, usecols=(0, 1, 3), unpack=True)
    
    # k, pk0, pk2, pk4 = pk_inst.pspec_old(x_c, y_c, z_c)

    # np.savetxt("./test/1.041halo.dat_LOS985_pow_old.gcat", np.array([k, pk0, pk2, pk4]).T)


    pk_inst.pspec_new(x_c, y_c, z_c, "/1.041halo.dat_LOS985_pypow_new.gcat")

def plot_aux(ax, file, usecols, color, label):
    k, pk0, pk2, pk4 = np.loadtxt(file, usecols=usecols, unpack=True)

    range_ = (k >=0.02) & (k <=0.5)
    k, pk0, pk2, pk4  = k[range_], pk0[range_], pk2[range_], pk4[range_]
    k_p = np.linspace(0.025, 0.495, 48)

    ax[0].plot(k_p, pk0, color=color, label=label)
    ax[1].plot(k_p, pk2, color=color, label=label)
    ax[2].plot(k_p, pk4, color=color, label=label)

    return k_p, pk0, pk2, pk4

def plot_aux_err(ax, file, color, label):
    k, pk0, pk2, pk4, std0, std2, std4 = np.loadtxt(file, usecols=(0, 1, 2, 3, 4, 5, 6), unpack=True)

    range_ = (k >=0.02) & (k <=0.5)
    k, pk0, pk2, pk4, std0, std2, std4 = k[range_], pk0[range_], pk2[range_], pk4[range_], std0[range_], std2[range_], std4[range_]
    k_p = np.linspace(0.025, 0.495, 48)

    ax[0].errorbar(k_p, pk0, yerr=std0, color=color, label=label)
    ax[1].errorbar(k_p, pk2, yerr=std2, color=color, label=label)
    ax[2].errorbar(k_p, pk4, yerr=std4, color=color, label=label)

    return std0* np.sqrt((0.5 ** 3) / 24) , std2* np.sqrt((0.5 ** 3) / 24) , std4* np.sqrt((0.5 ** 3) / 24)
    
def plot():
    import matplotlib.pyplot as pt

    fig, ax = pt.subplots(2, 3, figsize=(15, 10), sharex=True)

    k__, pkr0, pkr2, pkr4 = plot_aux(ax[0], "/global/project/projectdirs/desi/cosmosim/cov_mat_challenge/cubic_box/SLICS/PK_multipoles_V2/SLICS_Nmesh512_pkl_rsd/985.txt", (0,1,2,3), "orange", "SLICS_V2")
    std0, std2, std4 = plot_aux_err(ax[0], "/global/homes/a/avariu/desi_project_dir/FastPM_SLICS/scripts/avg_std_slics_pspec.dat", color="black", label="avg_slics_v2")

    kr, pko0, pko2, pko4 = plot_aux(ax[0], "./test/1.041halo.dat_LOS985_pow_old.gcat", (0, 1, 2, 3), "red", "pow_old")
    ax[1][0].plot(kr, (pkr0 - pko0) / std0, color="red")
    ax[1][1].plot(kr, (pkr2 - pko2) / std2, color="red")
    ax[1][2].plot(kr, (pkr4 - pko4) / std4, color="red")

    kp, pkp0, pkp2, pkp4 = plot_aux(ax[0], "./test/1.041halo.dat_LOS985.gcat_powspec_cheng.pspec", (0, 5, 6, 7), "cyan", "powspeccheng")
    ax[1][0].plot(kr, (pkr0 - pkp0) / std0, color="cyan", ls="-")
    ax[1][1].plot(kr, (pkr2 - pkp2) / std2, color="cyan", ls="-")
    ax[1][2].plot(kr, (pkr4 - pkp4) / std4, color="cyan", ls="-")


    kn, pkn0, pkn2, pkn4 = plot_aux(ax[0], "./test/1.041halo.dat_LOS985_pypow_new.gcat", (0, 5, 6, 7), "blue", "pypow_new")
    ax[1][0].plot(kr, (pkr0 - pkn0) / std0, color="blue", ls="--")
    ax[1][1].plot(kr, (pkr2 - pkn2) / std2, color="blue", ls="--")
    ax[1][2].plot(kr, (pkr4 - pkn4) / std4, color="blue", ls="--")

    
    kc, pkc0, pkc2, pkc4 = plot_aux(ax[0], "/global/project/projectdirs/desi/cosmosim/cov_mat_challenge/cubic_box/SLICS/PK_multipoles_V2_clean/SLICS_Nmesh512_pkl_rsd/985.txt", (0,1,2,3), "orange", "SLICS_V2_clean")
    ax[1][0].plot(kr, (pkr0 - pkc0) / std0, color="orange")
    ax[1][1].plot(kr, (pkr2 - pkc2) / std2, color="orange")
    ax[1][2].plot(kr, (pkr4 - pkc4) / std4, color="orange")

    ax[0][0].set_xscale("log")
    ax[0][1].set_xscale("log")
    ax[0][2].set_xscale("log")

    ax[1][1].set_ylim([-0.1, 0.1])
    ax[1][0].set_ylabel("$\\sigma$")

    ax[0][0].legend()
    fig.savefig("./test/powspec.png")

if __name__ == '__main__':
    # main()
    plot()
