import numpy as np
import matplotlib.pyplot as plt
import sys
from pylab import meshgrid

kf = open(sys.argv[1])
ifaces = np.array( kf.readline().split(), dtype=np.float ) 
q = np.array( kf.readline().split(), dtype=np.float ) 
q = np.append( q, (q[-1]+q[-2]) ) # for pcolor mesh
Y,X = meshgrid( ifaces, q )   

K = np.loadtxt(sys.argv[1], comments="#", skiprows=2)
nx, ny = np.shape(K)

fig = plt.figure( )
ax1 = fig.add_axes( [.10,.10,.35,.80] )
ax2 = fig.add_axes( [.5,.10,.35,.80], sharex=ax1, sharey=ax1 )
axcb = fig.add_axes( [.9,.10,.05,.80] )

# Real plot
ax1.pcolormesh(X, Y, K[0:nx//2].T , cmap="RdBu", vmin=-np.max(np.abs(K)), vmax=np.max(np.abs(K)) )
ax1.set_ylim( ifaces[-1], ifaces[0] )
ax1.set_xlim( q[-1], q[0] )
ax1.set_xscale("log", nonposx='clip')
#plt.colorbar()

# imaginary 
im = ax2.pcolormesh(X, Y, K[nx//2:].T , cmap="RdBu", vmin=-np.max(np.abs(K)), vmax=np.max(np.abs(K)) )
plt.setp(ax2.get_yticklabels(), visible=False)
plt.colorbar(im, axcb)

ax1.set_ylabel("Depth (m)")
ax1.set_xlabel("Pulse moment (A$\cdot$s)")
ax2.set_xlabel("Pulse moment (A$\cdot$s)")

plt.show()
exit()


plt.matshow(K[0:nx//2])
plt.colorbar()

plt.matshow(K[nx//2:])
plt.colorbar()

KA = np.abs(K[0:nx//2] + 1j*K[nx//2:])
plt.matshow(1e9*KA, cmap='viridis')
plt.colorbar()
plt.show()
