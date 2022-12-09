import numpy as np


#open text file in read mode
# text_file = open("output/out.para_1296_1000_0.5_64_48_m_DIF_cov_vdisp_pk.out", "r")
text_file = open("./out.test.out", "r")
 
#read whole file to a string
data = text_file.read()

substr = "chi2 = "
res = [i for i in range(len(data)) if data.startswith(substr, i)]
print(len(res))

substr = "in"

res_f = []
for i in range(len(res)):
    for k in range(len(data[res[i]: res[i] + 150])):
        if data[res[i]: res[i] + 150].startswith(substr, k):
            res_f.append(k)

print(len(res_f))

chi2_arra = np.zeros(len(res))

for i, k in enumerate(res):
    substr = data[k+7: res[i] + res_f[i]]
    if substr != "1e+101":
        chi2_arra[i] = float(substr)

min_indx_a = np.where(np.min(chi2_arra[chi2_arra>0]) == chi2_arra)
print(min_indx_a)

min_indx = min_indx_a[0][0]
print(min_indx, res[min_indx], data[res[min_indx] -8: res[min_indx] + res_f[min_indx]])

max_indx = np.argmax(chi2_arra)
print(max_indx, res[max_indx], data[res[max_indx] -8: res[max_indx] + res_f[max_indx]])
#close file

print(len(res), len(chi2_arra[chi2_arra > 0]))
text_file.close()