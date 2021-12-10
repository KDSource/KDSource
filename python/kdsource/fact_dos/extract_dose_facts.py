import numpy as np

filename = "log_dose_facts.txt"
file = open(filename, "w")
file.close()
file = open(filename, "a")

Es,fs = np.loadtxt("ARN_neutron").T
Es *= 1E-6
np.savetxt(file, [np.log(Es),np.log(fs)], fmt='%.2e', delimiter=',')
file.write("\n")

Es,fs = np.loadtxt("ARN_photon").T
Es *= 1E-6
np.savetxt(file, [np.log(Es),np.log(fs)], fmt='%.2e', delimiter=',')
file.write("\n")

Es,fs = np.loadtxt("ICRP_neutron").T
Es *= 1E-6
np.savetxt(file, [np.log(Es),np.log(fs)], fmt='%.2e', delimiter=',')
file.write("\n")

import matplotlib.pyplot as plt
log_Es = np.log(Es)
log_fs = np.log(fs)
plt.plot(log_Es,log_fs)
plt.show()
log_Es[0] = -1E3
log_fs[0] = -1E3
plt.plot(log_Es,log_fs)
plt.show()

Es,fs = np.loadtxt("ICRP_photon").T
Es *= 1E-6
np.savetxt(file, [np.log(Es),np.log(fs)], fmt='%.2e', delimiter=',')
file.write("\n")