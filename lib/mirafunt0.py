'''This module defines the analytic function, the varous numeric functions and the mean square differences between the values achieved for all. The numeric functions implement a finite time for decay.'''
from __future__ import division
import os.path
import numpy, random


N = 5# input('Enter Number of Cycles N (5): ')
GAMMA = 1# input('Enter value of GAMMA (1.0): ')
OMEGA_1 = 5#  input('Enter value of OMEGA_1(5.0): ')

################################### INITIAL PARAMETERS ###############################################
# DEFINE PRIMARY IMPORTANT CONSTANTS
'''N = 5 # how many full natural-oscialltion periods to model for
GAMMA = 1 # decay parameter in units of lambda
OMEGA_1 = 5.0 # Rabi frequancy in units of lambda'''

# CALCULATE MORE IMPORTANT CONSTANTS FROM ABOVE VALUES 
GAMMA_T = GAMMA/2
OMEGA = abs(OMEGA_1**2-((GAMMA_T-GAMMA)**2)/4)**0.5
k1 = GAMMA*GAMMA_T / (GAMMA*GAMMA_T + OMEGA_1**2)
k2 = OMEGA_1**2 / (GAMMA*GAMMA_T)
w1 = (GAMMA + GAMMA_T) / 2
k3 = (GAMMA + GAMMA_T) / (2*OMEGA)
T = (2*numpy.pi / OMEGA_1) * N # the full time the model will run for
#################################### END OF PARAMETERS ###############################################

class mainfunt0:
    #################################### SAVING FUNCTION #################################################
     
    @staticmethod
    def isdir(f):
        ''' Checks for directory and if needed makes it. Enter full path i.e. "C:\Users\owner" '''
        if not os.path.isdir(f):
            os.mkdir(f)

    #################################### SAVING FUNCTION #################################################
    @staticmethod
    def save(dirname,filename,data):
        savepath = dirname
        mainfunt0.isdir(savepath)
        completeName = os.path.join(savepath,filename)
        numpy.savetxt(completeName,data)
        print "saved " + filename
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    
    ############################## HERE STARTS THE ANALYTIC FUNCTION #####################################
    @staticmethod
    def W(t):
        '''returns the value of w which is defined in the source paper as w = p_aa - p_bb'''
        w = k1*(1+k2*numpy.exp(-w1*(t))*( numpy.cos(OMEGA*t) + k3*numpy.sin(OMEGA*t)))
        return w

    @staticmethod
    def anal_prob(t):
        '''returns the analytic probability of the particle being found in the excited state'''
        p_bb = 0.5*(1-mainfunt0.W(t))
        return p_bb
    
    @staticmethod
    def analpop(n):
        '''returns an array of how the prability of being in the excited state changes over the time period'''
        DT = 2*numpy.pi/(OMEGA_1*n)
        time = numpy.arange(0,T, DT)
        waveform = numpy.zeros(len(time))
        for i in range(len(time)):
            waveform[i] = mainfunt0.anal_prob(time[i])
        Wavename = "AnalyticSOL_GAMMA" +str(GAMMA) + "_OMEGA" + str(OMEGA_1) +"_N" + str(N) + "_n" + str(n)
        Timename = Wavename + "TIME"
        mainfunt0.save("Analytic Solutions", Wavename, waveform) 
        mainfunt0.save("Analytic Solutions", Timename, time) 
        return waveform, time
    
    ############################## HERE ENDS THE ANALYTIC FUNCTION #######################################






    ############################## HERE STARTS THE FIRST NUMERIC FUNCTION ######################################
    ############################### THIS USES THE SIMPLE SIN^2 MAJIGGER ########################################
    @staticmethod
    def f1(step, DT):
        '''returns the next co-ordinate in the cos plot'''
        theta = OMEGA_1/2*step*DT # this is how far through the period of 2PI it is
	rho = numpy.sin(theta)**2 # z is the probability of the qubit being in the excited state
	return rho
    
    @staticmethod
    def do_decay1(step, DT, stationarysol):
        '''randomises whether the qubit decays in this period'''
        if stationarysol[step]*GAMMA*DT > random.random():
            return 1
        else:
            return step + 1
            
    @staticmethod
    def get_next_position1(step, DT, stationarysol):
        '''returns the next postion in the wavepath of the single qubit along with it's new step number'''
        new_step = mainfunt0.do_decay1(step, DT, stationarysol)
        new_position = stationarysol[new_step]
        return (new_position, new_step)

    @staticmethod
    def single_waveform1(DT, stationarysol,number_of_steps):
        '''uss the previous functions to create a full wavepath for one qubit'''
        waveform = numpy.zeros((number_of_steps,2))
        for i in range(1,number_of_steps):
            waveform[i] = mainfunt0.get_next_position1(waveform[i-1,1], DT, stationarysol)
            # number_of_decays = len(waveform[:,1][waveform[:,0]==0]) - 1 # finds the number of decays in the time period. The -1 takes into account that it essentially starts in the decayed state
        return waveform[:,0]
    @staticmethod
    def sinpop(n,Q):
        '''uses the previous functions to create a population waveform of Q qubits with resolution of n data points in a full period'''
        #set up initial values 
        DT = 2*numpy.pi/(OMEGA_1*n)# time step in units of 1 / (rabi frequency)
        number_of_steps = n*N
        time = numpy.arange(0,T, DT)
        
        # data stores
        stationarysol = numpy.zeros(number_of_steps)
        wavepaths = [0]*Q
        average_wavepath = numpy.zeros(number_of_steps)
        
        # lets create a full undecaying f(t) for the time period of interest so we can access it later
        for i in range(number_of_steps):
                stationarysol[i] = mainfunt0.f1(i, DT)
    
    # creates Q individual wavepaths 
        for i in range(Q):
            wavepaths[i] = mainfunt0.single_waveform1(DT, stationarysol, number_of_steps)
    
    # sums up the wavepaths to create a population waveform
        for i in range(number_of_steps):
            total = 0 # this will add up all the combined probabilities of being in excited state
            for j in range(Q):
                total += wavepaths[j][i] 
            average_wavepath[i] = total/Q # this returns the average prob of being in excited state at this time step
            
        Wavename = "SINLSOL_GAMMA" +str(GAMMA) + "_OMEGA" + str(OMEGA_1) +"_N" + str(N) + "_n" + str(n)
        Timename = Wavename + "TIME"
        mainfunt0.save("Sin Solutions T0", Wavename, average_wavepath) 
        mainfunt0.save("Sin Solutions T0", Timename, time) 
        return average_wavepath, wavepaths[0], time
    ############################## HERE ENDS THE FIRST NUMERIC FUNCTION ########################################







    ############################## HERE STARTS THE SECOND NUMERIC FUNCTION ######################################
    ############################# THIS USES THE EULER METHOD OF INTEGRATION #####################################
    @staticmethod
    def do_decay2(ab, DT):
        '''randomises whether the qubit decays in this period'''
        pbb = (ab[1]*ab[1].conjugate()).real
        if pbb*GAMMA*DT > random.random():
            return 0 # decays
        else:
            return 1 # doesn't decay

    @staticmethod
    def take_step2((a,b),dt):
        '''takes a single simple integration step'''
        da = - (OMEGA_1/2 * b * dt)*1j
        db = - (OMEGA_1/2 * a * dt)*1j
    
        return (a+da, b+db)
   
    @staticmethod 
    def normalise2( (a,b ) ):
        '''Normailses the 2 relative probabilities fed in and returns the normalised values'''
        scale = abs(a)**2 + abs(b)**2
        a_norm = a / scale
        b_norm = b / scale
        return (a_norm, b_norm)

    @staticmethod
    def single_wavepath2((a0,b0), n):
        '''uses the previous functions to map the wavepath of a single qubit with random decay'''
        dt = 2*numpy.pi/(OMEGA_1*n)# time step in units of 1 / (rabi frequency)
        number_of_steps = n*N
        time = numpy.arange(0,T, dt)
        ab = numpy.zeros((number_of_steps,2)) + numpy.zeros((number_of_steps,2))*1j# a in [:;0] b in [:;1]
        ab[0] = (a0,b0)
        for i in range(1,number_of_steps): 
            if mainfunt0.do_decay2(ab[i-1],dt) == 0: 
                ab[i] = (1.0 + 0.0j, 0.0 + 0.0j)
        else:
                ab[i] = mainfunt0.normalise2(mainfunt0.take_step2(ab[i-1], dt))
        
        prob_ab = ab*(ab.conjugate()) 
   
        return prob_ab, time

    @staticmethod
    def single_wavepath_bitch2((a0,b0), number_of_steps,dt):
        '''uses the previous functions to map the wavepath of a single qubit with random decay'''
        ab = numpy.zeros((number_of_steps,2)) + numpy.zeros((number_of_steps,2))*1j# a in [:;0] b in [:;1]
        ab[0] = (a0,b0)
        for i in range(1,number_of_steps): 
            ab[i] = mainfunt0.normalise2(mainfunt0.take_step2(ab[i-1], dt))
            if mainfunt0.do_decay2(ab[i],dt) == 0: 
                ab[i] = (1.0 + 0.0j, 0.0 + 0.0j) # this bit ensure that it takes a period dt to decay
                
        
        prob_ab = ab*(ab.conjugate()) 
   
        return prob_ab[:,1] # returns only the b values

    @staticmethod
    def eulpop(n,Q):
        '''uses the previous functions to create a population waveform of Q qubits with resolution of n data points in a full period'''
        #set up initial values 
        dt = 2*numpy.pi/(OMEGA_1*n)# time step in units of 1 / (rabi frequency)
        number_of_steps = n*N
        time = numpy.arange(0,T, dt)
        (a0,b0) = (1.0 + 0.0j,0.0 + 0.0j)
        
        # data stores
        wavepaths = [0]*Q
        average_wavepath = numpy.zeros(number_of_steps)
    
        # creates Q individual wavepaths 
        for i in range(Q):
            wavepaths[i] = mainfunt0.single_wavepath_bitch2((a0,b0), number_of_steps, dt)
    
    # sums up the wavepaths to create a population waveform
        for i in range(number_of_steps):
            total = 0 # this will add up all the combined probabilities of being in excited state
            for j in range(Q):
                total += wavepaths[j][i] 
        average_wavepath[i] = total.real/Q # this returns the average prob of being in excited state at this time step
        
        Wavename = "EulerSOL_GAMMA" +str(GAMMA) + "_OMEGA" + str(OMEGA_1) +"_N" + str(N) + "_n" + str(n)
        Timename = Wavename + "TIME"
        mainfunt0.save("Euler Solutions T0", Wavename, average_wavepath) 
        mainfunt0.save("Euler Solutions T0", Timename, time) 
        return average_wavepath, time


    ############################## HERE ENDS THE SECOND NUMERIC FUNCTION ########################################


    ############################## HERE STARTS THE THIRD NUMERIC FUNCTION #######################################
    ########################## THIS USES THE RUNGE KUTTA METHOD OF INTEGRATION ##################################
    @staticmethod
    def f_rk(c):
        return - OMEGA_1/2 * c * 1.0j

    @staticmethod
    def take_step_rk((a,b),dt):
        '''takes a single simple integration step using my own runge-kutta method for the system'''
        k1a = mainfunt0.f_rk(b) ; k1b = mainfunt0.f_rk(a)
        k2a = mainfunt0.f_rk(b+(dt/2)*k1b) ; k2b = mainfunt0.f_rk(a+(dt/2)*k1a)
        k3a = mainfunt0.f_rk(b+(dt/2)*k2b) ; k3b = mainfunt0.f_rk(a+(dt/2)*k2a)
        k4a = mainfunt0.f_rk(b+dt*k3b) ; k4b = mainfunt0.f_rk(a+dt*k3a)
    
        da = dt/6 * (k1a + 2*k2a + 2*k3a + k4a)
        db = dt/6 * (k1b + 2*k2b + 2*k3b + k4b)
    
        return (a+da, b+db)
    
    @staticmethod
    def single_wavepath_rk((a0,b0), n):
        '''uses the previous functions to map the wavepath of a single qubit with random decay'''
        dt = 2*numpy.pi/(OMEGA_1*n)# time step in units of 1 / (rabi frequency)
        number_of_steps = n*N
        time = numpy.arange(0,T, dt)
        ab = numpy.zeros((number_of_steps,2)) + numpy.zeros((number_of_steps,2))*1j# a in [:;0] b in [:;1]
        ab[0] = (a0,b0)
        for i in range(1,number_of_steps): 
            if mainfunt0.do_decay2(ab[i-1],dt) == 0: 
                ab[i] = (1.0 + 0.0j, 0.0 + 0.0j)
            else:
                ab[i] = mainfunt0.normalise2(mainfunt0.take_step_rk(ab[i-1], dt))
        
        prob_ab = ab*(ab.conjugate()) 
   
        return prob_ab, time

    @staticmethod
    def single_wavepath_bitch_rk((a0,b0), number_of_steps,dt):
        '''copy of function above to be implemented in function below'''
        ab = numpy.zeros((number_of_steps,2)) + numpy.zeros((number_of_steps,2))*1j# a in [:;0] b in [:;1]
        ab[0] = (a0,b0)
        for i in range(1,number_of_steps):
            ab[i] = mainfunt0.normalise2(mainfunt0.take_step_rk(ab[i-1], dt))
            if mainfunt0.do_decay2(ab[i],dt) == 0: 
                ab[i] = (1.0 + 0.0j, 0.0 + 0.0j)  # this take no time to decay             
        prob_ab = ab*(ab.conjugate()) 
    
        return prob_ab[:,1] # returns only the b values

    @staticmethod
    def rkpop(n,Q):
        '''uses the previous functions to create a population waveform of Q qubits with resolution of n data points in a full period'''
        #set up initial values 
        dt = 2*numpy.pi/(OMEGA_1*n)# time step in units of 1 / (rabi frequency)
        number_of_steps = n*N
        time = numpy.arange(0,T, dt)
        (a0,b0) = (1.0 + 0.0j,0.0 + 0.0j)
            
        # data stores
        wavepaths = [0]*Q
        average_wavepath = numpy.zeros(number_of_steps)
        
        # creates Q individual wavepaths 
        for i in range(Q):
            wavepaths[i] = mainfunt0.single_wavepath_bitch_rk((a0,b0), number_of_steps, dt)
        
        # sums up the wavepaths to create a population waveform
        for i in range(number_of_steps):
            total = 0 # this will add up all the combined probabilities of being in excited state
            for j in range(Q):
                total += wavepaths[j][i] 
            average_wavepath[i] = total.real/Q # this returns the average prob of being in excited state at this time step
            
        Wavename = "T0_RungeKuttaSOL_GAMMA" +str(GAMMA) + "_OMEGA" + str(OMEGA_1) +"_N" + str(N) + "_n" + str(n) + "_Q" + str(Q)
        Timename = Wavename + "TIME"
        mainfunt0.save("Runge Kutta Solutions T0", Wavename, average_wavepath) 
        mainfunt0.save("Runge Kutta Solutions T0", Timename, time) 
        return average_wavepath, time
    ############################## HERE ENDS THE THIRDS NUMERIC FUNCTION ########################################
    
    
    
    
    
    
    
    ############################# ANY EXTRA FUNCTIONS GO HERE FOR NOW ####################################
    @staticmethod
    def msd(a,b):
        '''as a fraction'''
        A = numpy.array(a)
        B = numpy.array(b)
        m_s_d = sum(( ( (A[1:] - B[1:]) / A[1:] )**2 )) / len(A[1:])
        return m_s_d

    @staticmethod  
    def compare_Q_sin(n,Qs):
        analwave, analtime = mainfunt0.analpop(n)
        msd = numpy.zeros((len(Qs),2))
        msd[:,1] = Qs
        for i in range(len(Qs)):
            numericwave, singlewave, numerictime = mainfunt0.sinpop(n,Qs[i])
            msd[i,0] = mainfunt0.msd(analwave, numericwave)
        msdname = "SinMSDQs" + str(Qs) + "_GAMMA" +str(GAMMA) + "_OMEGA" + str(OMEGA_1) +"_N" + str(N)
        mainfunt0.save("Sin Solutions T0", msdname, msd) 
        return msd

    @staticmethod
    def compare_Q_eul(n,Qs):
        analwave, analtime = mainfunt0.analpop(n)
        msd = numpy.zeros((len(Qs),2))
        msd[:,1] = Qs
        for i in range(len(Qs)):
            numericwave, numerictime = mainfunt0.eulpop(n,Qs[i])
            msd[i,0] = mainfunt0.msd(analwave, numericwave)
        msdname = "EulerMSDQs" + str(Qs) + "_GAMMA" +str(GAMMA) + "_OMEGA" + str(OMEGA_1) +"_N" + str(N)
        mainfunt0.save("Euler Solutions T0", msdname, msd) 
        return msd

    @staticmethod    
    def compare_Q_rk(n,Qs):
        analwave, analtime = mainfunt0.analpop(n)
        msd = numpy.zeros((len(Qs),2))
        msd[:,1] = Qs
        for i in range(len(Qs)):
            numericwave, numerictime = mainfunt0.rkpop(n,Qs[i])
            msd[i,0] = mainfunt0.msd(analwave[0:len(numericwave)], numericwave)
        msdname = "RungeKuttaMSD_Qs" + str(Qs) + "_GAMMA" +str(GAMMA) + "_OMEGA" + str(OMEGA_1) +"_N" + str(N)
        mainfunt0.save("Runge Kutta Solutions T0", msdname, msd) 
        return msd

    @staticmethod    
    def compare_n_rk(ns,Q):
        msd = numpy.zeros((len(ns),2))
        msd[:,1] = ns
        for i in range(len(ns)):
            analwave, analtime = mainfunt0.analpop(ns[i])
            numericwave, numerictime = mainfunt0.rkpop(ns[i],Q)
            msd[i,0] = mainfunt0.msd(analwave[0:len(numericwave)], numericwave)
        msdname = "RungeKuttaMSD_ns" + str(ns) + "_GAMMA" +str(GAMMA) + "_OMEGA" + str(OMEGA_1) +"_N" + str(N)
        mainfunt0.save("Runge Kutta Solutions T0", msdname, msd) 
        return msd 
        
    ####################################### END OF EXTRA #################################################
    
    
if __name__ == "__main__" :
    n = 1000; QS = [10, 100,1000,10000] 
    mainfunt0.compare_Q_rk(n,QS)