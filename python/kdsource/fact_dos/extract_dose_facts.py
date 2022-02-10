import matplotlib.pyplot as plt

import numpy as np

filename = "log_dose_facts.txt"
file = open(filename, "w")
file.close()
file = open(filename, "a")

Es, fs = np.loadtxt("ARN_neutron").T
Es *= 1e-6
np.savetxt(file, [np.log(Es), np.log(fs)], fmt="%.2e", delimiter=",")
file.write("\n")

Es, fs = np.loadtxt("ARN_photon").T
Es *= 1e-6
np.savetxt(file, [np.log(Es), np.log(fs)], fmt="%.2e", delimiter=",")
file.write("\n")

Es, fs = np.loadtxt("ICRP_neutron").T
Es *= 1e-6
np.savetxt(file, [np.log(Es), np.log(fs)], fmt="%.2e", delimiter=",")
file.write("\n")

log_Es = np.log(Es)
log_fs = np.log(fs)
plt.plot(log_Es, log_fs)
plt.show()
log_Es[0] = -1e3
log_fs[0] = -1e3
plt.plot(log_Es, log_fs)
plt.show()

Es, fs = np.loadtxt("ICRP_photon").T
Es *= 1e-6
np.savetxt(file, [np.log(Es), np.log(fs)], fmt="%.2e", delimiter=",")
file.write("\n")
