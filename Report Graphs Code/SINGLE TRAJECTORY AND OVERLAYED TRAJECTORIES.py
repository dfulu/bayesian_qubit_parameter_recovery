from __future__ import division
import matplotlib.pyplot as pyplot
import numpy
import os.path

def load(dirname,filename):
        savepath = r"C:\Users\User\Documents\FIZZIX\3rd Year\Comp Project\Data" + dirname
        completeName = os.path.join(savepath,filename)
        data = numpy.loadtxt(completeName)
        return data



RKDIR = "\Runge Kutta Solutions"
RKMCPATHNAME = "RK_MC_population_[Q,n,GAMMA,RF,N] = [16384, 200, 1.0, 5.0, 25]"
RKPOP = load(RKDIR,RKMCPATHNAME)

ALLRKNAME = "ALLRKS_[Q,n,GAMMA,RF,N] = [1024, 200, 1.0, 5.0, 25]"
ALLRK = load(RKDIR,ALLRKNAME)


ANALDIR = "\Analytic Solutions"
ANALPOPNAME = "Analytic_population_[n,GAMMA,RF,N] = [200, 1.0, 5.0, 25]"
ANALPOP = load(ANALDIR, ANALPOPNAME)

ANALTIMENAME = "Analytic_population_[n,GAMMA,RF,N] = [200, 1.0, 5.0, 25]TIME"
ANALTIME = load(ANALDIR, ANALTIMENAME)

RKTIMENAME = "RKTIME_[n,N] = [200, 25]"
RKTIME = load(RKDIR, RKTIMENAME)

R0 =( (ANALPOP[1:len(RKPOP)] - RKPOP[1:])/ANALPOP[1:len(RKPOP)] )

MSDNAME = "RungeKuttaMSDQs_[n,GAMMA,RF,N] = [200, 1.0, 5.0, 25]__Qs = [1, 4, 16, 64, 256, 1024, 4096, 16384]"
x = load(RKDIR,MSDNAME)
MSD = x[:,0]
Qs = x[:,1]

ALPHA = 0.07
scale = 3
padding = 0.05
FS = 9 * scale
x0, y0 = 0.1, 0.85

def cm2inch(value):
    return value/2.54

FIGSIZE=(cm2inch(10*scale), cm2inch(6*scale))


pyplot.figure(figsize = FIGSIZE)
#PLOT THE WAVEFUNCTIONS
ax1 = pyplot.subplot2grid((2,1), (0,0))
pyplot.plot(RKTIME, ALLRK[1], color = 'red')
pyplot.ylim(-padding,1+padding)
ax1.text(x0, y0, r'(a)', fontsize=FS)
pyplot.ylabel(r'$b^{ 2} (t)$', fontsize=FS*1.4, labelpad=0*scale)
pyplot.yticks(numpy.arange(0, 1.5, 0.5), fontsize = FS)
pyplot.xticks(fontsize = FS)

ax2 = pyplot.subplot2grid((2,1), (1,0),sharex = ax1)
for i in range (100): 
    pyplot.plot(RKTIME, ALLRK[i], color = 'black', alpha = ALPHA,linewidth = 1.5)
pyplot.ylabel(r'$b^{ 2} (t)$', fontsize=FS*1.5, labelpad=0*scale)
pyplot.xlabel('Time [in units $1 / $'+r'$\Gamma$]', fontsize=FS, labelpad=-2*scale)
pyplot.ylim(-padding,1+padding)
ax2.text(x0, y0, r'(b)', fontsize=FS)
pyplot.yticks(numpy.arange(0, 1.5, 0.5), fontsize = FS)
pyplot.xticks(fontsize = FS)
pyplot.show()