import numpy as np
from matplotlib import pyplot as plt

q = 10000

ANALPROB = np.loadtxt("ANALPROB.txt")
ITFREQ = np.loadtxt("ITFREQ.txt")
CENTRES = np.loadtxt("CENTRES.txt")
RKPOP = np.loadtxt("RKPOP.txt")
ANALPOP = np.loadtxt("ANALPOP.txt")
R0 = np.loadtxt("R0.txt")

ANALTIME = np.loadtxt("ANALTIME.txt")
RKTIME = np.loadtxt("RKTIME.txt")
MODWAIT = np.loadtxt("MODWAIT.txt")


fig = plt.figure()
ax1 = plt.subplot2grid((2,3), (0,0), colspan = 2)
ax1.plot(ANALTIME, ANALPOP, color = 'black', label = 'Analytic solution')
ax1.plot(RKTIME, RKPOP, color = 'red', label = 'Numeric solution')
plt.legend(loc='best', fancybox=True, framealpha=0)
plt.ylabel('Population in state 1')
ax1.set_yticklabels(['','0.1','','0.3','', '0.5','', '0.7','', '0.9'])
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='on',         # ticks along the top edge are off
    labelbottom='off')

ax2 = plt.subplot2grid((4,3), (2,0), colspan = 2)
ax2.plot(RKTIME,R0, color = 'red')
plt.ylabel('Residual - %')
plt.xlabel('Time in units 1/'+r'$\Omega$' )
ax2.set_yticklabels(['','-10','','-6','', '-2','', '2','', '6'])



ax5 = plt.subplot2grid((4,3), (2, 2), sharey = ax2)
ax5.tick_params(axis='both', labelbottom='off', labelleft='off', left='off')
ax5.hist(R0, bins = 20,orientation='horizontal', color='red', normed = True,histtype='stepfilled')
ax5.text(0.14, -9, 'Residual\ndistribution', style='italic',)
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off')


plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
plt.show()
