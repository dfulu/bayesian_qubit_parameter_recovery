'''This module defines the analytic function, the varous numeric functions and the mean square differences between the values achieved for all. The numeric functions implement a finite time for decay.'''
from __future__ import division
import os.path
import numpy, random
import matplotlib.pyplot as pyplot  
import matplotlib.colors
import matplotlib.cm
import lib.isdir as isdir
GAMMA = 1 


def pp(t, RFsG): 
        '''Delivers the probability of the qubit decaying within time t for the given RF and the global GAMMA. Returns single value'''
        l = ((RFsG**2-(GAMMA/2)**2)**0.5)/2
        pt = GAMMA**3/(16*l**2) * ( (-2/GAMMA * numpy.sin(l*t) - 8*l/(GAMMA**2)*numpy.sin(l*t)*numpy.cos(l*t))*numpy.exp(-GAMMA*t/2) + 16*l**2/(GAMMA**3) * (1-numpy.exp(-GAMMA/2 * t)) )
        return pt
   
RFsG = numpy.exp(numpy.arange(0.01, 12,0.01))
TIMES= numpy.pi/RFsG

PPT = pp(TIMES, RFsG)

s = 1.4
scale = 3
padding = 0.05
FS = 9 * scale

def cm2inch(value):
    return value/2.54    

FIGSIZE=(cm2inch(10*scale), cm2inch(6*scale))
ax1 = pyplot.figure(figsize = FIGSIZE)
pyplot.subplots_adjust(
left  = 0.08,  # the left side of the subplots of the figure
right = 0.99,    # the right side of the subplots of the figure
)    

pyplot.loglog(RFsG, PPT, color = 'r')
#pyplot.xlabel(r'$\Omega_{0}$ / $\Gamma$', fontsize = FS)
#pyplot.ylabel('$q$', fontsize = FS*1.4)
pyplot.tick_params(which='major', reset=False, width = 2, length = 8)
pyplot.xticks(fontsize = FS*s)
pyplot.yticks(fontsize = FS*s)

pyplot.plot(numpy.asarray([-1, max(RFsG)]), numpy.asarray([0,0]) , color='black')
pyplot.show()