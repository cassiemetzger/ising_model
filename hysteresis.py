import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("H-data.csv", delimiter=",")

H = data[:, 0]
average_magnetization = data[:, 1]
analytical_phase_transition = 2.269

plt.title("H vs. Average Magnetization")
plt.xlabel("H")
plt.ylabel(r'|$\langle M \rangle$|')
plt.plot(H, average_magnetization, marker='x', linestyle='', label='Average Magnetization')
plt.savefig("H-magnetization.png")
plt.close()