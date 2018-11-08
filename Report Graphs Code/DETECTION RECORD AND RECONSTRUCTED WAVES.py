'''This module defines the analytic function, the varous numeric functions and the mean square differences between the values achieved for all. The numeric functions implement a finite time for decay.'''
from __future__ import division
import os.path
import numpy, random
import matplotlib.pyplot as pyplot  
import matplotlib.colors
import matplotlib.cm
import lib.isdir as isdir

N = 25# input('Enter Number of Cycles N (5): ')
GAMMA = 1.0# input('Enter value of GAMMA (1.0): ')
#OMEGA_1 = #  input('Enter value of OMEGA_1(5.0): ')
# CALCULATE MORE IMPORTANT CONSTANTS FROM ABOVE VALUES 
 # the full time the model will run for
#################################### END OF PARAMETERS ###############################################

class mainfunt0:
    #################################### SAVING FUNCTION #################################################
    @staticmethod
    def save(dirname,filename,data):
        savepath = r"C:\Users\User\Documents\FIZZIX\3rd Year\Comp Project\Data" + dirname
        isdir.isdir.isdir(savepath)
        completeName = os.path.join(savepath,filename)
        numpy.savetxt(completeName,data)
        print 'saved ' + filename +'\n'
    
    @staticmethod
    def load(dirname,filename):
        savepath = r"C:\Users\User\Documents\FIZZIX\3rd Year\Comp Project\Data" + dirname
        completeName = os.path.join(savepath,filename)
        data = numpy.loadtxt(completeName)
        return data

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    
    ############################## HERE STARTS THE ANALYTIC FUNCTION #####################################
    @staticmethod
    def W(t,k1,k2,k3,w1,OMEGA):
        '''returns the value of w which is defined in the source paper as w = p_aa - p_bb'''
        w = k1*(1+k2*numpy.exp(-w1*t)*( numpy.cos(OMEGA*t) + k3*numpy.sin(OMEGA*t)))
        return w

    @staticmethod
    def anal_prob(t,k1,k2,k3,w1,OMEGA):
        '''returns the analytic probability of the particle being found in the excited state'''
        p_bb = 0.5*(1-mainfunt0.W(t,k1,k2,k3,w1,OMEGA))
        return p_bb
    
    @staticmethod
    def analpop(n,RF):
        '''returns an array of how the prability of being in the excited state changes over the time period'''
        OMEGA_1 = RF
        GAMMA_T = GAMMA/2
        OMEGA = abs(OMEGA_1**2-((GAMMA_T-GAMMA)**2)/4)**0.5 
        k1 = GAMMA*GAMMA_T / (GAMMA*GAMMA_T + OMEGA_1**2)
        k2 = OMEGA_1**2 / (GAMMA*GAMMA_T) 
        w1 = (GAMMA + GAMMA_T) / 2
        k3 = (GAMMA + GAMMA_T) / (2*OMEGA)
        T = (1 / OMEGA_1) * N
        
        DT = 1/(OMEGA_1*n)
        time = numpy.arange(0,T, DT)
        waveform = numpy.zeros(len(time))
        for i in range(len(time)):
            waveform[i] = mainfunt0.anal_prob(time[i],k1,k2,k3,w1,OMEGA)
        Wavename = "Analytic_population_[n,GAMMA,RF,N] = " + str([n,GAMMA,RF,N])
        Timename = Wavename + "TIME"
        mainfunt0.save("\Analytic Solutions", Wavename, waveform) 
        mainfunt0.save("\Analytic Solutions", Timename, time) 
        return waveform, time
    
    ############################## HERE ENDS THE ANALYTIC FUNCTION #######################################





    ############################## HERE STARTS THE SECOND NUMERIC FUNCTION ######################################
    ############################# THIS USES THE EULER METHOD OF INTEGRATION #####################################

    @staticmethod
    def take_step_eul(RF,(a,b),dt):
        '''takes a single simple integration step'''
        da = - (RF/2 * b * dt)*1.0j
        db = - (RF/2 * a * dt)*1.0j - GAMMA * dt/2 * b
    
        return (a+da, b+db)
   
    @staticmethod 
    def normalise2( (a,b) ):
        '''Normailses the 2 relative probabilities fed in and returns the normalised values'''
        scale = abs(a)**2 + abs(b)**2
        a_norm = a / scale
        b_norm = b / scale
        return (a_norm, b_norm)

    @staticmethod
    def single_wavepath_bitchEUL(RF, (a0,b0), number_of_steps,dt):
        '''copy of function above to be implemented in function below'''
        JUMP = random.random();   P = [1,1]
        ab = numpy.zeros((number_of_steps,2)) + numpy.zeros((number_of_steps,2))*1j# a in [:;0] b in [:;1]
        ab[0] = (a0,b0)
        for i in range(1,number_of_steps):
            
            if (P[-1]<JUMP and P[0] > JUMP): 
                
                ddt = dt - ( (1-JUMP/P[0])/(GAMMA*abs(ab[i-1,1])**2) ) # find time in time step left after jump
                P[-1] = (1-ddt*GAMMA*abs(ab[i-1,1])**2) # correct the new P for after jump
                JUMP = random.random()
                ab[i] = mainfunt0.normalise2(mainfunt0.take_step_eul(RF, ab[0], ddt))
                P.append(P[-1]*(1-dt*GAMMA*abs(ab[i,1])**2)) 
                P = P[-2:] # Update P but keep previous value in 1st element of list
                               
            else: 
                ab[i] = mainfunt0.normalise2(mainfunt0.take_step_eul(RF, ab[i-1], dt)) 
                P.append(P[-1]*(1-dt*GAMMA*abs(ab[i,1])**2)); P = P[-2:] # Update P but keep previous value in 1st element of list
        
        prob_ab = abs(ab*(ab.conjugate()))
   
        return prob_ab[:,1] # returns only the b values

               

    @staticmethod
    def eulpop(n,Q,RF):
        '''uses the previous functions to create a population waveform of Q qubits with resolution of n data points in a full period'''
        #set up initial values 
        dt = 1/(RF*n)# time step in units of 1 / (rabi frequency)
        T = (1 / RF) * N
        number_of_steps = n*N
        time = numpy.arange(0,T, dt)
        (a0,b0) = (1.0 + 0.0j,0.0 + 0.0j)
        
        # data stores
        wavepaths = numpy.zeros((Q,number_of_steps))
        average_wavepath = numpy.zeros(number_of_steps)
    
        # creates Q individual wavepaths 
        for i in range(Q):
            wavepaths[i] = mainfunt0.single_wavepath_bitchEUL(RF, (a0,b0), number_of_steps, dt)
    
        # sums up the wavepaths to create a population waveform
        for i in range(number_of_steps):
             average_wavepath[i] = numpy.average(wavepaths[:,i])
        
        Wavename = "EulerSOL_[Q,n,GAMMA,RF,N] = " + str([Q,n,GAMMA,RF,N])
        Timename = Wavename + "TIME"
        mainfunt0.save("\Euler Solutions", Wavename, average_wavepath) 
        mainfunt0.save("\Euler Solutions", Timename, time) 
        return average_wavepath, time


    ############################## HERE ENDS THE SECOND NUMERIC FUNCTION ########################################




    ############################## HERE STARTS THE THIRD NUMERIC FUNCTION #######################################
    ########################## THIS USES THE RUNGE KUTTA METHOD OF INTEGRATION ##################################
        
    @staticmethod
    def fa_rk((a,b), RF):
        return - RF/2 * b * 1.0j
     
    @staticmethod   
    def fb_rk((a,b), RF):
        return - RF/2 * a * 1.0j - GAMMA/2 * b

    @staticmethod
    def take_step_rk((a,b),dt,RF):
        '''takes a single simple integration step using my own runge-kutta method for the system'''
        k1a = mainfunt0.fa_rk((a,b),RF) ; k1b = mainfunt0.fb_rk((a,b),RF)
        k2a = mainfunt0.fa_rk((a+(dt/2)*k1a,b+(dt/2)*k1b),RF) ; k2b = mainfunt0.fb_rk((a+(dt/2)*k1a,b+(dt/2)*k1b),RF)
        k3a = mainfunt0.fa_rk((a+(dt/2)*k2a,b+(dt/2)*k2b),RF) ; k3b = mainfunt0.fb_rk((a+(dt/2)*k2a,b+(dt/2)*k2b),RF)
        k4a = mainfunt0.fa_rk((a+dt*k3a,b+dt*k3b),RF) ; k4b = mainfunt0.fb_rk((a+dt*k3a,b+dt*k3b),RF)
    
        da = dt/6 * (k1a + 2*k2a + 2*k3a + k4a)
        db = dt/6 * (k1b + 2*k2b + 2*k3b + k4b)
    
        return (a+da, b+db)
   
    @staticmethod 
    def normaliseRK( (a,b ),dt):
        '''Normailses the 2 relative probabilities fed in and returns the normalised values'''
        scale = ( (abs(a)**2 + abs(b)**2))**0.5        
        a_norm = a / scale
        b_norm = b / scale
        return (a_norm, b_norm)

    @staticmethod
    def single_wavepath_rk(RF, n, N):
        '''generates single wavepath for frequency given'''
        dt = 1/(RF*n)# time step in units of 1 / (rabi frequency)
        number_of_steps = n*N
        time = numpy.arange(0,number_of_steps)*dt
        JUMP = random.random();   P = [1,1]
        ab = numpy.zeros((number_of_steps,2)) + numpy.zeros((number_of_steps,2))*1j# a in [:;0] b in [:;1]
        ab[0] = (1.0 + 0.0j,0.0 + 0.0j)
        decaytimes = [] # store the times where it has decayed
        for i in range(1,number_of_steps):
            
            if P[-1]<JUMP: 
                
                ddt = dt - ( (1-JUMP/P[0])/(GAMMA*abs(ab[i-1,1])**2) ) # find residual time in time step left after jump
                decaytimes.append(i*dt - ddt)
                P[-1] = (1-ddt*GAMMA*abs(ab[i-1,1])**2) # correct the new P for after jump
                JUMP = random.random()
                ab[i] = mainfunt0.normaliseRK(mainfunt0.take_step_rk(ab[0], ddt,RF),ddt)
                P.append(P[-1]*(1-dt*GAMMA*abs(ab[i,1])**2)); P = P[-2:] # Update P but keep previous value in 1st element of list
                                 
            else: 
                ab[i] = mainfunt0.normaliseRK(mainfunt0.take_step_rk(ab[i-1], dt,RF),dt) 
                P.append(P[-1]*(1-dt*GAMMA*abs(ab[i,1])**2)); P = P[-2:] # Update P but keep previous value in 1st element of list
               
                                            
        #generate returnable data
        DECT = numpy.asarray(decaytimes)
        prob_ab = abs(ab*(ab.conjugate()))
    
        return (prob_ab[:,1], time, DECT) # returns only the b values plus the time and the times of detection
        
   
    @staticmethod
    def wave_recon(RF, time, DECT):
        '''generates single wavepath from detection record'''
        dt = time[1] # time step in units of 1 / (rabi frequency)
        ab = numpy.zeros((len(time),2)) + numpy.zeros((len(time),2))*1j# a in [:;0] b in [:;1]
        ab[0] = (1.0 + 0.0j,0.0 + 0.0j)
        DNo = 0
        for i in range(1,len(time)):
                
            if  time[i] > DECT[DNo] and time[i-1]< DECT[DNo]:
                ddt = time[i] - DECT[DNo]                 
                ab[i] = mainfunt0.normaliseRK(mainfunt0.take_step_rk(ab[0], ddt,RF),ddt)
                if len(DECT) != DNo + 1:
                    DNo +=1     
            else: 
                ab[i] = mainfunt0.normaliseRK(mainfunt0.take_step_rk(ab[i-1], dt,RF),dt) 
        prob_ab = (ab*(ab.conjugate()) ).real
    
        return (prob_ab[:,1]) # returns only the b values
       
        

    @staticmethod
    def single_wavepath_bitch_rk((a0,b0), number_of_steps,dt,RF):
        '''copy of function above to be implemented in function below'''
        JUMP = random.random();   P = [1,1]
        ab = numpy.zeros((number_of_steps,2)) + numpy.zeros((number_of_steps,2))*1j# a in [:;0] b in [:;1]
        ab[0] = (a0,b0)
        decaytimes = [0] # store the times where it has decayed
        for i in range(1,number_of_steps):
            
            if (P[-1]<JUMP and P[0] > JUMP): 
                
                ddt = dt - ( (1-JUMP/P[0])/(GAMMA*abs(ab[i-1,1])**2) ) # find time in time step left after jump
                decaytimes.append(i*dt - ddt)
                P[-1] = (1-ddt*GAMMA*abs(ab[i-1,1])**2) # correct the new P for after jump
                JUMP = random.random()
                ab[i] = mainfunt0.normaliseRK(mainfunt0.take_step_rk(ab[0], ddt,RF),ddt)
                P.append(P[-1]*(1-dt*GAMMA*abs(ab[i,1])**2)) 
                P = P[-2:] # Update P but keep previous value in 1st element of list
                               
            else: 
                ab[i] = mainfunt0.normaliseRK(mainfunt0.take_step_rk(ab[i-1], dt,RF),dt) 
                P.append(P[-1]*(1-dt*GAMMA*abs(ab[i,1])**2)); P = P[-2:] # Update P but keep previous value in 1st element of list
               
                                            
        #generate returnable data
        DECT = numpy.asarray(decaytimes)
        waittimes = DECT[1:] - DECT[0:-1] # generates the wait times from the decay times
        prob_ab = abs(ab*(ab.conjugate()))
    
        return (prob_ab[:,1], waittimes) # returns only the b values


    @staticmethod
    def rkpop(n,Q,RF):
        '''uses the previous functions to create a population waveform of Q qubits with resolution of n data points in a full period'''
        #set up initial values 
        dt = 1/(RF*n)# time step in units of 1 / (rabi frequency)
        T = (1 / RF) * N
        number_of_steps = n*N
        
        (a0,b0) = (1.0 + 0.0j,0.0 + 0.0j)
        
        
        # data stores
        wavepaths = numpy.zeros((Q,number_of_steps))
        average_wavepath = numpy.zeros(number_of_steps)
        time = numpy.arange(0,T, dt)
        waittimes = numpy.zeros(1)
        
        # creates Q individual wavepaths 
        for i in range(Q):
            x = mainfunt0.single_wavepath_bitch_rk((a0,b0), number_of_steps, dt,RF)
            wavepaths[i] = x[0]
            waittimes = numpy.concatenate((waittimes,x[1]))
        
        # sums up the wavepaths to create a population waveform
        for i in range(number_of_steps):
            average_wavepath[i] = numpy.average(wavepaths[:,i])

            
        Wavename = "RK_MC_population_[Q,n,GAMMA,RF,N] = " + str([Q,n,GAMMA,RF,N])
        Timename = "RKTIME_[n,N] = " + str([n,N])
        Waitname = "Waittimes_[Q,n,GAMMA,RF,N] = " + str([Q,n,GAMMA,RF,N])
        Allrksname = "ALLRKS_[Q,n,GAMMA,RF,N] = " + str([Q,n,GAMMA,RF,N])
        mainfunt0.save("\Runge Kutta Solutions", Allrksname, wavepaths)
        mainfunt0.save("\Runge Kutta Solutions", Waitname, waittimes) 
        mainfunt0.save("\Runge Kutta Solutions", Wavename, average_wavepath) 
        mainfunt0.save("\Runge Kutta Solutions", Timename, time) 
        return average_wavepath, time, waittimes
    ############################## HERE ENDS THE THIRDS NUMERIC FUNCTION ########################################
    
 
    
    ######## FUNCTIONS TO FIND THE VALUE OMEGA GIVEN THE VALUE GAMMA AND A DATA SET OF DECAY TIMES ##############
    
    
    
    @staticmethod
    def prob_Rfreq(dn,RF):
        '''takes a wait time dn and returns the probability of measuring it for the value of rabi frequency RF given'''
        lamda = ((RF**2-(GAMMA/2)**2)**0.5)/2
        w = GAMMA*(RF/(2*lamda))**2*numpy.sin(lamda*dn)**2*numpy.exp(-GAMMA*dn/2)
        return w   
        
    @staticmethod
    def bayes(dn,RFs,pn):
        '''Uses bayes theorem to calculate the new probabilities given the prior probabilities pn of a range of rabi frequencies RFs'''
        PN = numpy.zeros(len(RFs)) # store the new, not normailsed probs
        for i in range(len(RFs)):
            PN[i] = mainfunt0.prob_Rfreq(dn,RFs[i])*pn[i]
        TotalProb = sum(PN)
        return PN/TotalProb
    
    @staticmethod
    def bayeteration(dns, RFs):
        ProbConv = numpy.zeros((len(RFs),len(dns)+1)) # generate store for prob data ProbConv[:,i] gives the probs after considering i wait times
        ProbConv[:,0] = 1/len(RFs) # this sets all frequencies equally likely before looking at data
        for i in range (len(dns)):
            ProbConv[:,i+1] = mainfunt0.bayes(dns[i],RFs,ProbConv[:,i]) # generates convergence of probabilities 
        return ProbConv
    
    

    
    
    ###################################### END OF BAYESIAN JUNK #################################################
    
    
   ############################## Distribution Junk can go in here #############################################
    
    @staticmethod
    def FreqCalc(Xmin, Xmax, data, bins, norm):
        '''this function generates either a normalised or non-normalised frequency distribution. 
        norm = 1 means it will be normailsed. norm = 0 means not normaised '''
        
        
        width = (Xmax-Xmin)/bins
        freq = numpy.zeros(bins)
        centres = numpy.arange(0.5*width, bins*width + 0.5*width, width)
        
        for i in range(bins):
            MIN = i*width; MAX = (i+1)*width
            freq[i] = len(data[data>MIN][data[data>MIN]<MAX]) # count the number of values in this interval
        if norm == 1:
            freq = freq/(sum(freq)*width)
            normed = 'yes'
        else: 
            normed = 'no'
        
        filename = "Numeric_wait_distribution[#times, bins, normed] = " + str([len(data),bins,normed])
        centrename = "Numeric_wait_centres[#times, bins] = " + str([len(data),bins])
        mainfunt0.save("\Runge Kutta Solutions", filename, freq)
        mainfunt0.save("\Runge Kutta Solutions", centrename, centres)    
        return (centres, freq)
    
    @staticmethod
    def itfreqcalc(Xmin, Xmax, data, bins, norm):
        width = (Xmax-Xmin)/bins
        freqsets = [0]*len(data)
        centres = numpy.arange(0.5*width, bins*width + 0.5*width, width)
        for i in range(len(data)):
            dataset = data[:i]
            freqsets[i] = numpy.zeros(bins)
            for j in range(bins):
                MIN = j*width; MAX = (j+1)*width
                freqsets[i][j] = len(dataset[dataset>MIN][dataset[dataset>MIN]<MAX]) # count the number of values in this interval
            if norm == 1:
                freqsets[i] = freqsets[i]/(sum(freqsets[i])*width)
        
        return (centres, freqsets)
    
    
    @staticmethod
    def analprob(RF, Xmin, Xmax, N):
        width = (Xmax-Xmin)/N
        prob = numpy.zeros(N)
        centres = numpy.arange(0.5*width, N*width + 0.5*width, width)
        for i in range(N): 
            prob[i] = mainfunt0.prob_Rfreq(centres[i],RF)
        return (centres, prob)

    @staticmethod
    def dectfreqgenerator(RF, n, N, q):
        dects = numpy.zeros(0)
        for i in range(q):
            (waste, time, minidec) = mainfunt0.single_wavepath_rk(RF, n, N)
            dects = numpy.concatenate((dects,minidec))
        
        Xmin = 0;
        Xmax = time[-1]
        bins = n*N
        data = dects
        dectfreq = mainfunt0.FreqCalc(Xmin, Xmax, data, bins, 0) # zero for not normalised
        
        return dectfreq #should return fequency values of decay along the evolution
    
            
        
    ################################## END OF DISTRIBUTION JUNK #################################################
    
    
    ############################# ANY EXTRA FUNCTIONS GO HERE FOR NOW ####################################
    @staticmethod
    def msd(a,b):
        '''as a fraction'''
        A = numpy.array(a)
        B = numpy.array(b)
        m_s_d = sum(( ( (A[1:] - B[1:]) / A[1:] )**2 )) / len(A[1:])
        return m_s_d

    @staticmethod    
    def compare_Q_rk(n,Qs,RF):
        analwave, analtime = mainfunt0.analpop(n,RF)
        msd = numpy.zeros((len(Qs),2))
        msd[:,1] = Qs
        for i in range(len(Qs)):
            numericwave, numerictime,waittimes = mainfunt0.rkpop(n,Qs[i],RF)
            msd[i,0] = mainfunt0.msd(analwave, numericwave)
        msdname = "RungeKuttaMSDQs_[n,GAMMA,RF,N] = " + str([n,GAMMA,RF,N]) + "__Qs = " + str(Qs)
        mainfunt0.save("\Runge Kutta Solutions", msdname, msd) 
        return msd

    @staticmethod    
    def compare_n_rk(ns,Q,RF):
        msd = numpy.zeros((len(ns),2))
        msd[:,1] = ns
        for i in range(len(ns)):
            analwave, analtime = mainfunt0.analpop(ns[i])
            numericwave, numerictime = mainfunt0.rkpop(ns[i],Q,RF)
            msd[i,0] = mainfunt0.msd(analwave, numericwave)
        msdname = "RungeKuttaMSDns" + str(ns) + "_GAMMA" +str(GAMMA) + "_OMEGA" + str(RF) +"_N"
        mainfunt0.save("\Runge Kutta Solutions T0", msdname, msd) 
        return msd 
        
    ####################################### END OF EXTRA #################################################


    ##################################### PLOTTING FUNCTIONS #############################################
    @staticmethod
    def bayescolour(RFs,ITPROBS, n):
        '''this function plots the bayesian iteration method in a colour plot. It will plot the first n values of the iteration probs'''
        tick = 1 # can reset the RF tick frequency
        RFrange = numpy.ptp(RFs) #calculates the range of the rabi frequency values
        sca = n/RFrange
        pyplot.figure() 
        pyplot.imshow(ITPROBS[:,:n], origin = 'lower', cmap=pyplot.get_cmap('seismic'), extent = (0,n, min(RFs)*sca,max(RFs)*sca)) 
        pyplot.yticks(numpy.arange(min(RFs),max(RFs),tick)*sca, numpy.arange(min(RFs),max(RFs),tick));
        pyplot.show()
    
    @staticmethod
    def singlewavegraphs(RFs,n,N):
        #first generate the data
        wavepaths = [0] * len(RFs)
        for i in range(len(RFs)):
            if i == 0: 
                (wavepaths[i], time, DECT) = mainfunt0.single_wavepath_rk(RFs[i], n, N) 
            else:
                wavepaths[i] = mainfunt0.wave_recon(RFs[i], time, DECT)+i
            
        colour = ['red','blue','purple']
        pyplot.figure()
        for i in range(len(RFs)):
            pyplot.plot(time,wavepaths[i], color = colour[i])
        pyplot.yticks([0],[''])
        pyplot.show()
        
            
    
    
    ######################################## PLOTTING END ################################################
    
scale = 3
FS = 9 * scale*0.8
x0, y0 = 0.1, 0.7
pady = 0.05

def cm2inch(value):
    return value/2.54

FIGSIZE=(cm2inch(10*scale), cm2inch(6*scale))
#CREATE DATA AND PLOT
n = 200;  RF0 = 5.0; N = 50
RF = [RF0, RF0*2,RF0/2]
AX = [0,0,0,0]
WAVE =  [0,0,0]
WAVE[0], TIME, DECT = mainfunt0.single_wavepath_rk(RF[0], n, N)

COLOUR = ['red', 'blue','green']
STYLE = ['-','-','-']
WIDTH = 0.4*scale
LETTER=['a','b','c','d']
for i in range(2):
    WAVE[i+1] = mainfunt0.wave_recon(RF[i+1], TIME, DECT)



pyplot.figure()
ax1 = pyplot.subplot2grid((4,1), (0,0))
ax1.vlines(DECT,0, 1,linewidth = WIDTH)
ax1.tick_params(axis='both', labelbottom='off', labelleft='off', left='off', bottom = 'off')
ax1.text(x0, y0, '(' + LETTER[0] +')', fontsize=FS, bbox=dict(facecolor='white', edgecolor='none', pad=0, alpha = 0.8))

for i in range(3): 
    ax = pyplot.subplot2grid((4,1), (i+1,0), sharex = ax1)
    ax.plot(TIME,WAVE[i], color = COLOUR[i], linestyle = STYLE[i], linewidth = WIDTH)
    ax.tick_params(axis='both', labelbottom='off', labelleft='off')
    pyplot.ylim(0-pady,1+pady)
    ax.text(x0, y0, '(' + LETTER[i+1] +')', fontsize=FS, bbox=dict(facecolor='white', edgecolor='none', pad=0, alpha = 0.8))
ax.tick_params(axis='both', labelbottom='on')
pyplot.xlabel('Time [in units $1 / $'+r'$\Gamma$]', fontsize=FS, labelpad=0*scale)
pyplot.xticks(fontsize = FS)
pyplot.show() 