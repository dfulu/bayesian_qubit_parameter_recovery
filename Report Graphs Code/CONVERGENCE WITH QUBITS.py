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

ANALDIR = "\Analytic Solutions"
ANALPOPNAME = "Analytic_population_[n,GAMMA,RF,N] = [200, 1.0, 5.0, 25]"
ANALPOP = load(ANALDIR, ANALPOPNAME)

ANALTIMENAME = "Analytic_population_[n,GAMMA,RF,N] = [200, 1.0, 5.0, 25]TIME"
ANALTIME = load(ANALDIR, ANALTIMENAME)

RKTIMENAME = "RKTIME_[n,N] = [200, 25]"
RKTIME = load(RKDIR, RKTIMENAME)

R0 =( (ANALPOP[1:len(RKPOP)] - RKPOP[1:]) )#/ANALPOP[1:len(RKPOP)] )

MSDNAME = "RungeKuttaMSDQs_[n,GAMMA,RF,N] = [200, 1.0, 5.0, 25]__Qs = [1, 4, 16, 64, 256, 1024, 4096, 16384]"
x = load(RKDIR,MSDNAME)
MSD = x[:,0]
Qs = x[:,1]



scale = 3
padding = 0.05
FS = 9 * scale

def cm2inch(value):
    return value/2.54

FIGSIZE=(cm2inch(10*scale), cm2inch(5*scale))




pyplot.figure(figsize = FIGSIZE)

'''
#PLOT THE WAVEFUNCTIONS

pyplot.plot( RKTIME[0::STEP], RKPOP[0::STEP], color='blue', linestyle = '', marker = '.',markersize = 10)
pyplot.plot( ANALTIME, ANALPOP, color='#FF3399', label='analytic')
STEP = 75

pyplot.ylabel(r'$b^{ 2} (t)$', fontsize=FS, labelpad=0*scale)
pyplot.ylim(-padding + min(R0)*10 - 0.1, max(ANALPOP) + padding)
pyplot.xlim(-0.2, 5.2)
#pyplot.yticks(numpy.arange(0, 1, 0.2), fontsize = FS)
pyplot.xticks(fontsize = FS)

#PLOT THE RESIDUALS
pyplot.plot( RKTIME[1:], R0*10 - 0.1, color='blue')
pyplot.xlabel('Time [in units $1 / $'+r'$\Gamma$]', fontsize=FS, labelpad=0*scale)
pyplot.xticks(fontsize = FS)
pyplot.yticks(numpy.arange(0,1,0.2), fontsize = FS)

'''
#PLOT THE CONVERGENCE

ax4 = pyplot.subplot2grid((1,1), (0,0))
pyplot.plot(Qs, MSD, color='blue', linestyle = '', marker = '.',markersize = 20)
ax4.set_yscale('log');ax4.set_xscale('log')
pyplot.ylabel('RMSD', fontsize=FS, style = 'italic')
pyplot.xlabel('Number of Trajectories', fontsize=FS, style = 'italic')
pyplot.yticks(fontsize = FS)
pyplot.xticks(fontsize = FS)
#pyplot.legend(loc='best', fancybox=True, framealpha=0)
#pyplot.ylabel("Normalised MSD")


pyplot.show()