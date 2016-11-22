import numpy as np
import matplotlib.pyplot as plt
import sys

K = np.loadtxt(sys.argv[1], comments="#")
nx, ny = np.shape(K)

plt.matshow(K[0:nx//2])
plt.colorbar()

plt.matshow(K[nx//2:])
plt.colorbar()

KA = np.abs(K[0:nx//2] + 1j*K[nx//2:])
plt.matshow(1e9*KA, cmap='viridis')
plt.colorbar()
plt.show()
