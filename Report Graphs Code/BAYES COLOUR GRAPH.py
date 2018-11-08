from __future__ import division
import os.path, numpy
import matplotlib.pyplot as plt  
import matplotlib.colors
import matplotlib.cm
import lib.isdir as isdir
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
bayesname = "LIKELIHOOD UNDER UPDATE RULE FOR Waittimes_[Q,n,GAMMA,RF0,N,w] = [20000, 400, 1.0, 5.0, 75, 10.0]"
dfcentrename = "Numeric_wait_centres[#times, bins] = [373092, 2000] for Waittimes_[Q,n,GAMMA,RF0,N,w] = [20000, 400, 1.0, 5.0, 75, 10.0]"
BAYES = load(detuningdir,bayesname)
propdecaydir = r"\decay time distribution"
rfname = "RF0s_[n,GAMMA,w,N] = [400, 1.0, 10, 75]"
RF0s = load(propdecaydir,rfname)

PEAK = numpy.zeros(len(BAYES[1,:]-1))

for i in range(1,len(BAYES[1,:])): 
    index = numpy.where(BAYES[:,i] == max(BAYES[:,i]))[0][0]
    PEAK[i-1] = RF0s[index]
'''
def plot(n, BAYES,PEAK,RF0s):
    
    scale = 3
    padding = 0.05
    FS = 9 * scale

    def cm2inch(value):
        return value/2.54    
    

    FIGSIZE=(cm2inch(10*scale), cm2inch(6*scale))
    fig,im = plt.subplots(figsize = FIGSIZE)
    plt.subplots_adjust(
    left  = 0.08,  # the left side of the subplots of the figure
    right = 0.99,    # the right side of the subplots of the figure
    )    
    
    RFSMAX = max(RF0s); RFSMIN = min(RF0s)
    tick = 1 # can reset the RF tick frequency
    RFrange = numpy.ptp(RF0s) #calculates the range of the rabi frequency values
    sca = n/RFrange/3

    ax1 = plt.subplot2grid((1,1), (0,0))
    plt.ylabel(r'$\Omega_{0}$ [in units of $\Gamma$]', fontsize = FS)
    plt.xlabel('Number of decay times considered', fontsize = FS)
    im = ax1.imshow(BAYES[:,:n], origin = 'lower', cmap=plt.get_cmap('nipy_spectral'), extent = (0,n, RFSMIN*sca,RFSMAX*sca), norm=matplotlib.colors.LogNorm()) # 'Paired' 
    plt.yticks(numpy.arange(1,10,tick)*sca, numpy.arange(1,10,tick), fontsize = FS);
    plt.xticks(fontsize = FS);
    plt.ylim(0.05*sca,10*sca)
    plt.xlim(0,n)    
    plt.plot(numpy.arange(1,len(PEAK[:n])+1,1),PEAK[:n]*sca, color = 'black', lw = 2)
    ax1.tick_params(which='major', reset=False, width = 2, length = 8)
    
    
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("top", size="8%", pad=0.15)
    
    cb = plt.colorbar(im, cax=cax, orientation='horizontal',ticks=[10**(-60),10**(-50), 10**(-40),10**(-30), 10**(-20), 10**(-10), 1])#, ticks=[10**(-150), 10**(-120),10**(-90), 10**(-60), 10**(-30), 1])
    
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')    
    
    cb.ax.tick_params(labelsize=FS) 
    cb.set_label('P('+r'$\Omega$'+')', labelpad=10, y=0.5, fontsize = FS)
    cb.ax.tick_params(which='major', reset=False, width = 2, length = 8)
    #cbar = fig.colorbar(cax, ticks=[10**(-150), 10**(-100), 10**(-50), 1])
    
    #fig.colorbar(im, ax=ax1, )
    cb.ax.invert_xaxis()
    
    
    

    plt.show()
