import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
PROCO = np.loadtxt("PROCO")

ANALPROB = np.loadtxt("ANALPROB.txt")
ITFREQ = np.loadtxt("ITFREQ.txt")
RFS = np.loadtxt("RFS.txt")
CENTRES = np.loadtxt("CENTRES.txt")
RFSMAX = max(RFS)
RFSMIN = min(RFS)

'''i = 10; icol = 'pink'
j = 100; jcol = 'yellow' '''
n = 400



# First set up the figure, the axis, and the plot element we want to animate
fig,im = plt.subplots()


tick = 1 # can reset the RF tick frequency
RFrange = np.ptp(RFS) #calculates the range of the rabi frequency values
sca = n/RFrange/3

ax1 = plt.subplot2grid((1,1), (0,0))
plt.ylabel(r'$\Omega$')
plt.xlabel('Number of decay times considered')
im = ax1.imshow(PROCO[:,:n], origin = 'lower', cmap=plt.get_cmap('Paired'), extent = (0,n, RFSMIN*sca,RFSMAX*sca), norm=matplotlib.colors.LogNorm()) 
plt.yticks(np.arange(RFSMIN,RFSMAX,tick)*sca, np.arange(RFSMIN,RFSMAX,tick));
#ax1.set_ylim([RFSMIN, RFSMAX])
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.1)

cb = plt.colorbar(im, cax=cax, ticks=[10**(-150), 10**(-120),10**(-90), 10**(-60), 10**(-30), 1])
cb.set_label('P('+r'$\Omega$'+')', labelpad=0, y=0.5)
#cbar = fig.colorbar(cax, ticks=[10**(-150), 10**(-100), 10**(-50), 1])
#cbar.ax.set_yticklabels(['10**(-150)', '10**(-100)', '10**(-50)', '1'])
#fig.colorbar(im, ax=ax1, )
'''
plt.vlines(i,RFSMIN*sca,RFSMAX*sca, lw = 1, color = icol)
plt.vlines(j,RFSMIN*sca,RFSMAX*sca, lw = 1, color = jcol)

ax2 = plt.subplot2grid((1,3), (0,2), colspan = 1, sharey = ax1)
plt.yticks(np.arange(RFSMIN,RFSMAX,tick)*sca, np.arange(RFSMIN,RFSMAX,tick));

plt.plot(PROCO[:,i],(RFS-0.1)*sca, color = icol)
plt.plot(PROCO[:,j],(RFS-0.1)*sca, color = jcol)
'''
plt.show()
