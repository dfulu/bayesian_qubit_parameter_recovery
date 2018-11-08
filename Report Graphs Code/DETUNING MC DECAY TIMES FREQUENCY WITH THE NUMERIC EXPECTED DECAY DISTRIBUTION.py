from __future__ import division
import os.path
import numpy, random
import matplotlib.pyplot as pyplot  
import matplotlib.colors
import matplotlib.cm
import lib.isdir as isdir

def save(dirname,filename,data):
        savepath = r"C:\Users\User\Documents\FIZZIX\3rd Year\Comp Project\Data" + dirname
        isdir.isdir.isdir(savepath)
        completeName = os.path.join(savepath,filename)
        numpy.savetxt(completeName,data)
        print 'saved ' + filename +'\n'
    
def load(dirname,filename):
        savepath = r"C:\Users\User\Documents\FIZZIX\3rd Year\Comp Project\Data" + dirname
        completeName = os.path.join(savepath,filename)
        data = numpy.loadtxt(completeName)
        return data
'''
detuningdir = r"\Runge Kutta Solutions detuning"
decayfreqname = "Numeric_wait_distribution[#times, bins, normed] = [373092, 2000, 'yes'] for Waittimes_[Q,n,GAMMA,RF0,N,w] = [20000, 400, 1.0, 5.0, 75, 10.0]"
dfcentrename = "Numeric_wait_centres[#times, bins] = [373092, 2000] for Waittimes_[Q,n,GAMMA,RF0,N,w] = [20000, 400, 1.0, 5.0, 75, 10.0]"
propdecaydir = r"\decay time distribution"
propdecayname = "Decay_time_distribution[n,GAMMA,w,N] = [400, 1.0, 10, 75]"
rfname = "RF0s_[n,GAMMA,w,N] = [400, 1.0, 10, 75]"
pdtimename = "TIME[n,N] = [400, 75]"

DECAYFREQ = load(detuningdir,decayfreqname)
DFCENTRES = load(detuningdir,dfcentrename)
PROPDECAY = load(propdecaydir,propdecayname)
PDTIME = load(propdecaydir,pdtimename)
RF0s = load(propdecaydir, rfname)
'''


def plot(DECAYFREQ,DFCENTRES,PROPDECAY,PDTIME,RF0s):
    Rindex = numpy.where(RF0s == 5.0)[0][0]
    Tindex = numpy.where(PDTIME == 10.0)[0][0]
    T2index = numpy.where(DFCENTRES > 10.0)[0][0]
    
    scale = 3
    padding = 0.05
    FS = 9 * scale
    
    PROPDECAY = PROPDECAY[Rindex,:Tindex]
    PDTIME = PDTIME[:Tindex]
    
    
    
    
    def cm2inch(value):
        return value/2.54
    
    FIGSIZE=(cm2inch(10*scale), cm2inch(6*scale))




    pyplot.figure(figsize = FIGSIZE)
    pyplot.subplots_adjust(
    left  = 0.17,  # the left side of the subplots of the figure
    right = 0.97,    # the right side of the subplots of the figure
    bottom = 0.14,   # the bottom of the subplots of the figure
    top = 0.96,    # the top of the subplots of the figure
    wspace = 0.33,   # the amount of width reserved for blank space between subplots
    hspace = 0.33   # the amount of height reserved for white space between subplots
    )
    
    xmin = -0.2; xmax = 10.2
    
    #PLOT THE WAVEFUNCTIONS
    STEP = 3
    pyplot.plot( DFCENTRES[:T2index], DECAYFREQ[:T2index], color='#f59f02', linestyle = '', marker = '.',markersize = 10)#[0::STEP]
    pyplot.plot( PDTIME, PROPDECAY, color='blue', label='analytic')
    pyplot.hlines(0, -5, 10, colors='gray', linestyles='solid')
    pyplot.ylim(-0.05, 0.65)
    pyplot.xlim(xmin, xmax)
    #pyplot.yticks(numpy.arange(0, 1, 0.2), fontsize = FS)
    pyplot.xticks(fontsize = FS)
    
    #PLOT THE RESIDUALS
    pyplot.xlabel('Time [in units $1 / $'+r'$\Gamma$]', fontsize=FS, labelpad=0*scale)
    pyplot.xticks(fontsize = FS)
    pyplot.yticks(fontsize = FS)
    
    pyplot.ylabel(r'$\omega (t)$', fontsize=FS*1.3, labelpad = 6)
    
    pyplot.show()