import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("temperature-data.csv", delimiter=",")

energy = data[:, 0]
average_magnetization = data[:, 1]
average_susceptibility = data[:, 2]
average_heat_capacity = data[:, 3]
analytical_phase_transition = 2.269

plt.title("KT vs. Average Magnetization (J = 0.01)")
plt.xlabel("KT")
plt.ylabel(r'|$\langle M \rangle$|')
plt.plot(energy, average_magnetization, marker='x', linestyle='', label='Average Magnetization')
plt.axvline(x=analytical_phase_transition, color='red', linestyle='--')
plt.savefig("KT-magnetization-test9.png")
plt.close()

def percent_error(actual, theoretical):
    return (np.abs(actual - theoretical))/theoretical * 100

plt.title("KT vs. Average Susceptibility (J = 1.5)")
plt.xlabel("KT")
plt.ylabel(r'$\chi$')
plt.plot(energy, average_susceptibility, marker='x', linestyle='', label='Average Susceptibility')
plt.axvline(x=analytical_phase_transition, color='red', linestyle='--', label = 'Theoretical value = 2.269 J ')
idx = np.argmax(average_susceptibility)
maxx = energy[idx]
plt.axvline(x= maxx, color='green', linestyle='--', label = f'Actual value = {maxx} J')
plt.savefig("KT-susceptibility-test9.png")
plt.close()

print('Percent error for susceptibility: ', percent_error(maxx, analytical_phase_transition))

plt.title("KT vs. Average Heat Capacity (J = 1.5)")
plt.xlabel("KT")
plt.ylabel(r"$C_v$")
plt.plot(energy, average_heat_capacity, marker='x', linestyle='', label='Average Heat Capacity')
plt.axvline(x=analytical_phase_transition, color='red', linestyle='--', label = 'Theoretical value = 2.269 J ')
idx = np.argmax(average_heat_capacity)
maxx = energy[idx]
plt.axvline(x= maxx, color='green', linestyle='--', label = f'Actual value = {maxx} J')
plt.savefig("KT-heat-capacity-test9.png")
plt.close()
print('Percent error for heat capacity: ', percent_error(maxx, analytical_phase_transition))