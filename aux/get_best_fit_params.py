import numpy as np


def get_best_fit(filename=None, usecols=None):
    list_ = np.loadtxt(filename, usecols=usecols, unpack=True)

    idx = np.argmin(list_[0])

    text_ = ""
    for el_ in list_:
        text_ += f"{el_[idx]}, "


    print(text_)

path = "/global/homes/a/avariu/phd/hodsimulations/output/FastPM4SLICS_hod_results/"


# get_best_fit(filename=path+"/pk_diff_kmax_05/1296/MN_logMcut_11.6_13.6_sigma_logM_0.01_4.01_logM1_9.0_14.0_k_0.0_20.0_alpha_0.0_1.3_vdisp_0.7_1.5_.txt", usecols=(1,2,3,4,5,6,7))
# get_best_fit(filename=path+"/pk_diff_kmax_04/1296/MN_logMcut_11.6_13.6_sigma_logM_0.01_4.01_logM1_9.0_14.0_k_0.0_20.0_alpha_0.0_1.3_vdisp_0.7_1.5_.txt", usecols=(1,2,3,4,5,6,7))
# get_best_fit(filename=path+"/pk_diff_kmax_05_wospikes/1296/MN_logMcut_11.6_13.6_sigma_logM_0.01_4.01_logM1_9.0_14.0_k_0.0_20.0_alpha_0.0_1.3_vdisp_0.7_1.5_.txt", usecols=(1,2,3,4,5,6,7))
get_best_fit(filename=path+"/pk_diff_kmax_04_wospikes/1296/MN_logMcut_11.6_13.6_sigma_logM_0.01_4.01_logM1_9.0_14.0_k_0.0_20.0_alpha_0.0_1.3_vdisp_0.7_1.5_.txt", usecols=(1,2,3,4,5,6,7))

