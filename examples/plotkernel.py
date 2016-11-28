import numpy as np
import matplotlib.pyplot as plt
import sys
from pylab import meshgrid
from matplotlib.colors import LightSource
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import LogLocator

kf = open(sys.argv[1])
ifaces = np.array( kf.readline().split(), dtype=np.float ) 
q = np.array( kf.readline().split(), dtype=np.float ) 
q = np.append( q, (q[-1]+q[-2]) ) # for pcolor mesh
Y,X = meshgrid( ifaces, q )   

K = 1e9*np.loadtxt(sys.argv[1], comments="#", skiprows=2)
nx, ny = np.shape(K)

fig = plt.figure( )
ax1 =  fig.add_axes( [.100,.125,.335,.775] )
ax2 =  fig.add_axes( [.465,.125,.335,.775] , sharex=ax1, sharey=ax1 )
axcb = fig.add_axes( [.835,.125,.040,.775] )

ccmap = "coolwarm" # RdBu

# Real plot
ax1.pcolormesh(X, Y, K[0:nx//2].T, cmap=ccmap, vmin=-np.max(np.abs(K)), vmax=np.max(np.abs(K)))
ax1.set_ylim( ifaces[-1], ifaces[0] )
ax1.set_xlim( q[-1], q[0] )
ax1.set_xscale("log", nonposx='clip')
ax1.minorticks_off()
ax1.set_title(r"$\mathrm{Re} \left( \mathcal{K}_N^{1D} \right)$")

ax1.xaxis.set_major_locator(LogLocator(base = 2.0))
ax1.get_xaxis().set_major_formatter(ScalarFormatter())

# imaginary 
im = ax2.pcolormesh(X, Y, K[nx//2:].T, cmap=ccmap, vmin=-np.max(np.abs(K)), vmax=np.max(np.abs(K)))
plt.setp(ax2.get_yticklabels(), visible=False)
cb = plt.colorbar(im, axcb)
cb.set_label(r"$\overline{V}^{\left(0\right)}_N$  (nV)")

ax1.set_ylabel("Depth (m)")
ax1.set_xlabel(r"Pulse moment (A$\cdot$s)")
ax2.set_xlabel(r"Pulse moment (A$\cdot$s)")
ax2.set_title(r"$\mathrm{Im} \left( \mathcal{K}_N^{1D} \right)$")

plt.suptitle("exp 0")

plt.show()
exit()
