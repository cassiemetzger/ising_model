import numpy as np
import matplotlib.pyplot as plt

data2 = np.loadtxt("lattice-2KbT.csv", delimiter=",")
plt.title("Lattice at 2KbT")
plt.imshow(data2, cmap="viridis")
plt.colorbar()
plt.savefig("lattice2.png")
plt.close()

data3 = np.loadtxt("lattice-3KbT.csv", delimiter=",")
plt.title("Lattice at 3KbT")
plt.imshow(data3, cmap="viridis")
plt.colorbar()
plt.savefig("lattice3.png")
plt.close()

