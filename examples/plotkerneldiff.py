import numpy as np
import matplotlib.pyplot as plt
import sys

K1 = np.loadtxt(sys.argv[1], comments="#")
K2 = np.loadtxt(sys.argv[2], comments="#")
nx, ny = np.shape(K1)

plt.matshow( K1[0:nx//2] - K2[0:nx//2] , cmap='inferno')
plt.colorbar()

plt.matshow( K1[nx//2:] - K2[nx//2] , cmap='inferno')
plt.colorbar()

KAc =   K1[0:nx//2]+1j*K1[nx//2:]
KA2c =  K2[0:nx//2]+1j*K2[nx//2:]
#              (K1[0:nx//2]+1j*K1[nx//2:]) )
KD = np.abs(KAc) - np.abs(KA2c)

plt.matshow(KD, cmap='RdBu', vmin=-np.max(np.abs(KD)), vmax=np.max(np.abs(KD)) )
plt.colorbar()
plt.show()
