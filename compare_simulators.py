# Compare variables
import matplotlib.pyplot as plt
import numpy as np


data_es = np.genfromtxt("variables_neuron_simulator.dat")
data_ag = np.genfromtxt("eglif_CA1PC_1_multimeter-09-0.dat")



plt.figure()

plt.subplot(311)
plt.plot(data_es[:,0],data_es[:,1])
plt.plot(data_ag[:,1],data_ag[:,2])
plt.title("Vm")
plt.ylabel("Voltage [mV]")


plt.subplot(312)
plt.plot(data_es[:,0],data_es[:,2])
plt.plot(data_ag[:,1],data_ag[:,5])
plt.title("Iadap")
plt.ylabel("Current [pA]")

plt.subplot(313)
plt.plot(data_es[:,0],data_es[:,3])
plt.plot(data_ag[:,1],data_ag[:,4])
plt.title("Idep")
plt.xlabel("Time [ms]")
plt.ylabel("Current [pA]")
plt.legend(["python","pyNEST"])
plt.show()
