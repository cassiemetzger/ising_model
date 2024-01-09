import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("timestep-data-2.csv", delimiter=",")

time_steps = data[:, 0]
energy = data[:, 1]
magnetization = data[:, 2]

plt.title("Energy Time Data at 2KbT")
plt.xlabel("Timestep")
plt.ylabel("Energy")
plt.plot(time_steps, energy, label='Energy')
plt.tight_layout()
plt.savefig("energy-time-2.png")
plt.close()

plt.title("Magnetization Time Data at 2KbT")
plt.xlabel("Timestep")
plt.ylabel("Magnetization")
plt.plot(time_steps, magnetization, label='Magnetization')
plt.tight_layout()
plt.savefig("magnetization-time-2.png")
plt.close()

data = np.loadtxt("timestep-data-3.csv", delimiter=",")
time_steps = data[:, 0]
energy = data[:, 1]
magnetization = data[:, 2]

plt.title("Energy Time Data at 3KbT")
plt.xlabel("Timestep")
plt.ylabel("Energy")
plt.plot(time_steps, energy, label='Energy')
plt.tight_layout()
plt.savefig("energy-time-3.png")
plt.close()

plt.title("Magnetization Time Data at 3KbT")
plt.xlabel("Timestep")
plt.ylabel("Magnetization")
plt.plot(time_steps, magnetization, label='Magnetization')
plt.tight_layout()
plt.savefig("magnetization-time-3.png")
