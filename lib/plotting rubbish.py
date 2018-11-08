from __future__ import division
import matplotlib.pyplot as pyplot
import numpy

NW0 = numpy.loadtxt("RK_DT0_N1000_Q5000")
NWDT = numpy.loadtxt("RK_DTN_N1000_Q5000")
NUMTIME = numpy.loadtxt("RK_DTN_N1000_Q5000_TIME")

ANALWAVE = numpy.loadtxt("ANALWAVE_N1000")
ANALTIME = numpy.loadtxt("ANALWAVE_N1000_TIME")

R0 = ANALWAVE - NW0
RDT = ANALWAVE - NWDT


pyplot.figure()

ax1 = pyplot.subplot2grid((2,3), (0,0), colspan=3)
pyplot.plot(ANALTIME, ANALWAVE, color = "black", label = "analytic")
pyplot.plot(  color='red', label='$sin^2$')
pyplot.plot(NUMTIME,NW0,  color='red', label='instant decay')
pyplot.plot(NUMTIME,NWDT,  color='purple', label='DT decay Time')
pyplot.legend(loc='best', fancybox=True, framealpha=0)

#PLOT THE RESIDUALS
ax2 = pyplot.subplot2grid((2,3), (1,0), colspan=3)
pyplot.plot( ANALTIME, R0, color='red')
pyplot.plot( ANALTIME, RDT, color='purple')
pyplot.ylabel("Residuals - %")
#pyplot.xlim(0,6.2)


pyplot.show()