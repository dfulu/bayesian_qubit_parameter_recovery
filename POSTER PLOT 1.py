from __future__ import division
import matplotlib.pyplot as pyplot
import numpy

RKPOP = numpy.loadtxt("RKPOP.txt")
ANALPOP = numpy.loadtxt("ANALPOP.txt")
R0 = numpy.loadtxt("R0.txt")

ANALTIME = numpy.loadtxt("ANALTIME.txt")
RKTIME = numpy.loadtxt("RKTIME.txt")




MODWAIT = numpy.loadtxt("MODWAIT.txt")
ANALPROB = numpy.loadtxt("ANALPROB.txt")
CENTRES = numpy.loadtxt("CENTRES.txt")

RESIDUAL_1 = numpy.loadtxt("RESIDUAL_1n100q8000data.txt")
RESIDUAL_2 = numpy.loadtxt("RESIDUAL_2n100q8000data.txt") 
RESIDUAL_3 = numpy.loadtxt("RESIDUAL_3n100q8000data.txt") 
MSD1 = numpy.loadtxt("MSD1data.txt")
MSD2 = numpy.loadtxt("MSD2data.txt")
MSD3 = numpy.loadtxt("MSD3data.txt")
Qs = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]

pyplot.figure()
pyplot.suptitle("Evalutation of Monte Carlo Wavefunction", fontsize=20, style='italic')
#PLOT THE WAVEFUNCTIONS
ax1 = pyplot.subplot2grid((3,3), (0,0), colspan=2)
pyplot.plot( TIME, ANALWAVE, color='grey', label='analytic')
#pyplot.plot( TIME, NUMERICWAVE_1, color='red', label='numeric_1')
#pyplot.plot( TIME, NUMERICWAVE_2, color='pink', label='numeric_2')
pyplot.plot( TIME, NUMERICWAVE_3, color='purple', label='numeric_3')
pyplot.plot( RKDT0TIME, NUMERICWAVE_RK, color='black', label='numeric_dt0')
pyplot.ylabel("Excited Population")
pyplot.text(2, 0.2, 'Comparison for 8000 qubits', style='italic',)
pyplot.xlim(0,6.2)


#PLOT THE RESIDUALS
ax2 = pyplot.subplot2grid((3,3), (1,0), colspan=2)
pyplot.plot( TIME[1:], RESIDUAL_1, color='red')
pyplot.plot( TIME[1:], RESIDUAL_2, color='pink')
pyplot.plot( TIME[1:], RESIDUAL_3, color='purple')
pyplot.ylabel("Residuals - %")
pyplot.xlim(0,6.2)


#PLOT THE RESIDUAL HISTOGRAM
ax3 = pyplot.subplot2grid((3,3), (1, 2), sharey = ax2)
pyplot.tick_params(axis='both', labelbottom='off', labelleft='off', left='off')
pyplot.hist(RESIDUAL_3, bins = 20,orientation='horizontal', color='purple', normed = True,histtype='step')
pyplot.hist(RESIDUAL_2, bins = 20,orientation='horizontal', color='pink', normed = True,histtype='step')
pyplot.hist(RESIDUAL_1, bins = 20,orientation='horizontal', color='red', normed = True, histtype = 'step')
pyplot.text(0.14, -9, 'Residual\nFrequency', style='italic',)

#PLOT THE CONVERGENCE
ax4 = pyplot.subplot2grid((3,3), (2,2),)
pyplot.semilogx(Qs, MSD1, color='red', label = "MSD1")
pyplot.semilogx(Qs, MSD2, color='pink', label = "MSD2")
pyplot.semilogx(Qs, MSD3, color='purple', label = "MSD3")
pyplot.xlabel('Numer of Qubits')
pyplot.text(800, 0.003, 'Mean Square\nDifference', style='italic',)
#pyplot.legend(loc='best', fancybox=True, framealpha=0)
#pyplot.ylabel("Normalised MSD")


# PLOT SINGLE WAVEFORM
ax5 = pyplot.subplot2grid((3,3), (2,0), colspan = 2)
pyplot.plot( TIME, SINGLEWAVE_1, color='red')
pyplot.ylabel("Single Qubit Path")
pyplot.xlabel('Time in units 1/Rabi Period')
pyplot.xlim(0,6.2)



#PLOT LEGEND
ax6 = pyplot.subplot2grid((3,3), (0,2), frameon = False)
pyplot.plot(0,0,  color='red', label='$sin^2$')
pyplot.plot(0,0,  color='pink', label='$Euler$')
pyplot.plot(0,0,  color='purple', label='$Runge Kutta$')
pyplot.plot(0,0, color='grey', label='$Analytic$')
pyplot.legend(loc='center', fancybox=True, framealpha=0)
pyplot.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')

pyplot.show()