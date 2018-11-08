'''This module defines the analytic function, the varous numeric functions and the mean square differences between the values achieved for all. The numeric functions implement a finite time for decay.'''
from __future__ import division
import os.path
import numpy, random
import matplotlib.pyplot as pyplot  
import matplotlib.colors
import  matplotlib.cm

#import lib.isdir as isdir

N = 50# input('Enter Number of Cycles N (5): ')
GAMMA = 1.0# input('Enter value of GAMMA (1.0): ')
#OMEGA_1 = #  input('Enter value of OMEGA_1(5.0): ')

################################### INITIAL PARAMETERS ###############################################
# DEFINE PRIMARY IMPORTANT CONSTANTS
'''N = 5 # how many full natural-oscialltion periods to model for
GAMMA = 1 # decay parameter in units of lambda
OMEGA_1 = 5.0 # Rabi frequancy in units of lambda'''

# CALCULATE MORE IMPORTANT CONSTANTS FROM ABOVE VALUES 
 # the full time the model will run for
#################################### END OF PARAMETERS ###############################################

class mainfunt0:
    #################################### SAVING FUNCTION #################################################
    @staticmethod
    def save(dirname,filename,data):
        #savepath = r"C:\Users\owner\Documents\FIZZIX\3rd Year\Comp Project\Data" + dirname
        #isdir.isdir.isdir(savepath)
        #completeName = os.path.join(savepath,filename)
        #numpy.savetxt(completeName,data)
        print ''
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
        Wavename = "AnalyticSOL_GAMMA" +str(GAMMA) + "_RF" + str(OMEGA_1) +"_N" + str(N) + "_n" + str(n)
        Timename = Wavename + "TIME"
        mainfunt0.save("\Analytic Solutions", Wavename, waveform) 
        mainfunt0.save("\Analytic Solutions", Timename, time) 
        return waveform, time
    
    ############################## HERE ENDS THE ANALYTIC FUNCTION #######################################


    ############################## HERE STARTS THE THIRD NUMERIC FUNCTION #######################################
    ########################## THIS USES THE RUNGE KUTTA METHOD OF INTEGRATION ##################################
    @staticmethod
    def f_rk(c, RF):
        return - RF/2 * c * 1.0j

    @staticmethod
    def take_step_rk((a,b),dt,RF):
        '''takes a single simple integration step using my own runge-kutta method for the system'''
        k1a = mainfunt0.f_rk(b,RF) ; k1b = mainfunt0.f_rk(a,RF)
        k2a = mainfunt0.f_rk(b+(dt/2)*k1b,RF) ; k2b = mainfunt0.f_rk(a+(dt/2)*k1a,RF)
        k3a = mainfunt0.f_rk(b+(dt/2)*k2b,RF) ; k3b = mainfunt0.f_rk(a+(dt/2)*k2a,RF)
        k4a = mainfunt0.f_rk(b+dt*k3b,RF) ; k4b = mainfunt0.f_rk(a+dt*k3a,RF)
    
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
        prob_ab = (ab*(ab.conjugate())).real
    
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
        prob_ab = ab*(ab.conjugate()) 
    
        return (prob_ab[:,1], waittimes) # returns only the b values


    @staticmethod
    def rkpop(n,Q,RF):
        '''uses the previous functions to create a population waveform of Q qubits with resolution of n data points in a full period'''
        #set up initial values 
        dt = 1/(RF*n)# time step in units of 1 / (rabi frequency)
        T = (1 / RF) * N
        number_of_steps = n*N
        time = numpy.arange(0,T, dt)
        (a0,b0) = (1.0 + 0.0j,0.0 + 0.0j)
        waittimes = numpy.zeros(1)
        
        # data stores
        wavepaths = [0]*Q
        average_wavepath = numpy.zeros(number_of_steps)
        
        # creates Q individual wavepaths 
        for i in range(Q):
            x = mainfunt0.single_wavepath_bitch_rk((a0,b0), number_of_steps, dt,RF)
            wavepaths[i] = x[0]
            waittimes = numpy.concatenate((waittimes,x[1]))
        
        # sums up the wavepaths to create a population waveform
        for i in range(number_of_steps):
            total = 0 # this will add up all the combined probabilities of being in excited state
            for j in range(Q):
                total += wavepaths[j][i] 
            average_wavepath[i] = total.real/Q # this returns the average prob of being in excited state at this time step
            
        Wavename = "RungeKuttaSOL_GAMMA" +str(GAMMA) + "_OMEGA" + str(RF) +"_N" + str(N) + "_n" + str(n)
        Timename = Wavename + "TIME"
        mainfunt0.save("\Runge Kutta Solutions T0", Wavename, average_wavepath) 
        mainfunt0.save("\Runge Kutta Solutions T0", Timename, time) 
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
        analwave, analtime = mainfunt0.analpop(n)
        msd = numpy.zeros((len(Qs),2))
        msd[:,1] = Qs
        for i in range(len(Qs)):
            numericwave, numerictime = mainfunt0.rkpop(n,Qs[i],RF)
            msd[i,0] = mainfunt0.msd(analwave, numericwave)
        msdname = "RungeKuttaMSDQs" + str(Qs) + "_GAMMA" +str(GAMMA) + "_OMEGA" + str(RF)
        mainfunt0.save("\Runge Kutta Solutions T0", msdname, msd) 
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


'''savefig('demo.png', transparent=True)'''

 

n = 40; q =5000  ; RF = 5.0
ANAL, ANALTIME = mainfunt0.analpop(n,RF)
RK, TIME, WAIT = mainfunt0.rkpop(n,q,RF)
R0 =( (ANAL[1:len(RK)] - RK[1:])/ANAL[1:len(RK)] ) *100

pyplot.figure()
pyplot.suptitle("Evalutation of Monte Carlo Wavefunction", fontsize=20, style='italic')

#PLOT THE WAVEFUNCTIONS
ax1 = pyplot.subplot2grid((2,3), (0,0), colspan=3)
pyplot.plot( ANALTIME, ANAL, color='grey', label='Analytic')
pyplot.plot( TIME[:len(RK)], RK, color='purple', label='RK instant decay')
pyplot.ylabel("Excited Population")
#pyplot.xlim(0,6.2)

#PLOT THE RESIDUALS
ax2 = pyplot.subplot2grid((2,3), (1,0), colspan=3, sharex = ax1)
pyplot.plot( TIME[1:], R0, color='purple')
pyplot.ylabel("Residuals - %")
pyplot.xlabel("Time")
pyplot.xlim(0,6.2)

'''
pyplot.figure()
pyplot.suptitle("Evalutation of Monte Carlo Wavefunction", fontsize=20, style='italic')

#PLOT THE WAVEFUNCTIONS
ax1 = pyplot.subplot2grid((3,3), (0,0), colspan=2)
pyplot.plot( ANALTIME, ANAL, color='grey', label='Analytic')
pyplot.plot( TIME[:len(RK)], RK, color='purple', label='RK instant decay')
pyplot.ylabel("Excited Population")
#pyplot.xlim(0,6.2)


#PLOT THE RESIDUALS
ax2 = pyplot.subplot2grid((3,3), (1,0), colspan=2, sharex = ax1)
pyplot.plot( TIME[1:], R0, color='purple')
pyplot.ylabel("Residuals - %")
pyplot.xlabel("Time")
pyplot.xlim(0,6.2)


#PLOT THE RESIDUAL HISTOGRAM
ax3 = pyplot.subplot2grid((3,3), (1, 2), sharey = ax2)
pyplot.tick_params(axis='both', labelbottom='off', labelleft='off', left='off')
pyplot.hist(R0, bins = 20,orientation='horizontal', color='purple', normed = True,histtype='stepfilled')
pyplot.title("Residual Frequency", fontsize=12, style='italic')


# PLOT WAITITMES
ax4 = pyplot.subplot2grid((3,3), (2, 0), colspan = 2)
pyplot.tick_params(axis='both', labelbottom='on', labelleft='off', left='on')
pyplot.hist(WAIT[WAIT>0][WAIT[WAIT>0]<5], bins = 80,orientation='vertical', color='purple', normed = True,histtype='stepfilled')
pyplot.ylabel("Distribution")
pyplot.xlabel("Wait Time")
pyplot.text(3, 0.4, 'Wait Time\nFrequency', style='italic',)

#STATE n, q and OMEGA/GAMMA ON GRAPH
ax5 = pyplot.subplot2grid((3,3), (0,2), frameon = False)
ax5.text(0.5, 0.5, 'n = '+str(n) + '\nq = ' + str(q) + '\n' + r'$\Omega/\gamma$ = ' + str(OMEGA_1/GAMMA), fontsize = 12, style='italic',horizontalalignment='center',verticalalignment='center',transform=ax5.transAxes)
#pyplot.text(0, 0.5, 'n = '+str(n) + '\nq = ' + str(q) + '\n' + r'$\Omega/\gamma$ = ' + str(OMEGA_1/GAMMA), style='italic')
pyplot.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
ax5.set_xlim([0, 1]); ax5.set_ylim([0, 1])

ax6 = pyplot.subplot2grid((3,3), (2,2), frameon = False)
pyplot.text(0.5, 0.5, 'The Time throughout\nis in units of ' + r'$2\pi/\Omega$', fontsize = 12, style='italic',horizontalalignment='center',verticalalignment='center',transform=ax6.transAxes )
pyplot.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
ax6.set_xlim([0, 1]); ax6.set_ylim([0, 1])
'''
pyplot.show()
