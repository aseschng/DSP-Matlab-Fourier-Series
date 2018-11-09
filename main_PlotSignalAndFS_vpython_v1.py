# -*- coding: utf-8 -*-
"""
This module is to generate the Fourier Series coefficients of 4 different types
of signals.
@author: Chng Eng Siong
"""
import numpy as np
import matplotlib.pyplot as plt

def fn_getVariousSignals_FS_Coeff (typeSignal, numK, N, myFundamentalFreq, optionCell):
    """ The variables are:
    % para 1= 'DIY', 'Square', 'Saw', 'Triangle'
    % Para 2= number of sinusoids to use (includes DC) 
    % Para 3= number of samples in 1 sec (sampling freq)
    % Para 4= fundamental freqeuncy of the periodic signal to construct
    % Para 5= optionCell for various plots, e.g, DutyCycle for 'Square',
    %          'Ascending' == {1,0} for Saw
    """
    # y(t) = 10 +3cos(w_0 t) + 3cos(w_0 t + 0) + 5cos(2w_0+phi/6) + 4sin(3w_0 t-pi/2)
    if (typeSignal == 'DIY'):
        mySigFS = [(10,0*myFundamentalFreq,0), (3,myFundamentalFreq,0), (5,2*myFundamentalFreq,np.pi/6), (4,3*myFundamentalFreq,-np.pi/2)]
        return mySigFS



    # defining the coefficients of a square wave
    # a[n] = (2*A/(n*pi))*sin(n*pi/2)
    #% where n = 1,2,3,4...
    if (typeSignal == 'Square'):
       A = 1; DutyCycle = 0.5;
       if (optionCell[0] == 'DutyCycle'):
            DutyCycle   = optionCell[1]
            
      # The dc term
       mySigFS =  []
       mySigFS = mySigFS+[(A*DutyCycle, 0,0)]
       for k in np.arange(1,numK):
            tmpA   = (2*A/(k*np.pi))*np.sin(k*np.pi*DutyCycle)
            tmpF   = k*myFundamentalFreq
            tmpPhi  = 0
    
            if tmpA < 0:
                tmpA    = abs(tmpA);
                tmpPhi  = -np.pi;
                if (abs(tmpA)<1e-6):
                    tmpPhi = 0;

            mySigFS =  mySigFS+[(tmpA,tmpF,tmpPhi)]
    return mySigFS



    if (typeSignal == 'Triangle'):
        A = 1
        mySigFS  = []
        mySigFS  = mySigFS+[(A, 0, 0)]
        for k in np.arange(1,numK):
            if (k%2 == 1):
                tmpA    = ((8*A)/(np.pi*np.pi*k*k));
            else:
                tmpA    = 0                
                
            tmpF    =  k*myFundamentalFreq 
            tmpPhi  =  0
            mySigFS =  mySigFS+[(tmpA,tmpF,tmpPhi)]

        return mySigFS   # for Triangle Waves
    


    if (typeSignal == 'Saw'):
        A = 1
        mySigFS  = []
        mySigFS  = mySigFS+[(A/2, 0, 0)]
    
        for k in np.arange(1,numK):
            tmpA    = (2*A)/(k*2*np.pi)
            tmpF    = k*myFundamentalFreq 
            tmpPhi  =  +np.pi/2
            mySigFS =  mySigFS+[(tmpA,tmpF,tmpPhi)]

        return mySigFS   # for Saw Waves
    




def  fn_genTimeSignalFrom_FSCoeff(mySigFS, Fs):
    N  = Fs   # we generate 1 second worth of data
    w0 = 2*np.pi/N
    n  = np.arange(0,N)  # creating N samples, starting idx = 0
    t   = n/Fs
    K   = len(mySigFS)
    my_y = np.zeros((K,N))  #rows = each one signal from FS coeff with numSample=N
    my_yT = np.zeros(N)

    for idx,myOneFSCoeff in enumerate(mySigFS):
      my_y[idx,:] = myOneFSCoeff[0]*np.cos(w0*myOneFSCoeff[1]*n+myOneFSCoeff[2]);
      my_yT = my_yT+my_y[idx,:]
  
    return t,my_yT,my_y





def fn_PostProcessPhase(YT):

    magYT = abs(YT)/len(YT)
    maxYT = max(magYT);
    processedPhaseYT = np.angle(YT)*(magYT>(1/1000)*maxYT);
    # making the phase of the +ve frequency near pi to become -pi
    # making the phase of the -ve frequency near pi to become +pi

    halfN = np.ceil(len(YT)/2)
    A = processedPhaseYT;
    Arange = np.arange(0,len(YT));

    rangeLeftA = (abs((abs(A)-np.pi))<1e-6) & (Arange<halfN);
    A[rangeLeftA] = -np.pi;
    rangeRightA = (abs(abs(A)-np.pi)<1e-6) & (Arange>=halfN);
    A[rangeRightA] = +np.pi;
    
    rangeA = (abs(A-0)<1e-6);
    A[rangeA] = 0;
    processedPhaseYT = A;

    return(processedPhaseYT)    




"""
Lets generate the test signal's fourier coefficients
"""
Fs = 8000;   N= Fs;
#typeSignal = 'DIY'; numK = 4;  N = 8000; myFundamentalFreq = 10; optionCell = ()
typeSignal = 'Square'; numK = 35;   myFundamentalFreq = 10; optionCell = ('DutyCycle',0.1)
#typeSignal = 'Saw'; numK = 15;   myFundamentalFreq = 10
#typeSignal = 'Triangle'; numK = 15;   myFundamentalFreq = 10; optionCell = {}

tstSigFS      = fn_getVariousSignals_FS_Coeff (typeSignal, numK, N, myFundamentalFreq, optionCell)
print(tstSigFS)


K = len(tstSigFS);
highestFreq       = max( [x[1] for x in tstSigFS])
highestFreq       = min(Fs, highestFreq+2)



"""
Lets  generate the signals, plot each signal, its sum
"""
t, my_yT, my_y = fn_genTimeSignalFrom_FSCoeff(tstSigFS, Fs)
t10, my_yT10, my_y10 = fn_genTimeSignalFrom_FSCoeff(tstSigFS, Fs*10)

# Lets plot 3 cycles
numPeriodToPlot = 2
fundamentalFreq = min( [x[1] for x in tstSigFS[1:len(tstSigFS)]])
numsampleToPlot = int(np.ceil(numPeriodToPlot*(1/fundamentalFreq)*Fs));


fig, axs = plt.subplots(1,1)
for idx,myOneFSCoeff in enumerate(tstSigFS):
     myA = myOneFSCoeff[0]
     myFreq = myOneFSCoeff[1]
     myPhase = myOneFSCoeff[2]
     
     if (myPhase >= 0):
         tmpLegend =  "y{0:1d}(t) = {1:4.1e} cos(2\u03C0{2:3d}t + {3:3.2f} )".format(idx, myA,myFreq,myPhase)
     else:    
         tmpLegend =  "y{0:1d}(t) = {1:4.1e} cos(2\u03C0{2:3d}t  {3:3.2f} )".format(idx, myA,myFreq,myPhase)
  
     axs.plot(t[0:numsampleToPlot],my_y[idx,0:numsampleToPlot],label=tmpLegend)
         
axs.grid()
axs.set_xlabel('time(sec)')
axs.set_ylabel('y(t)')
plt.legend()
plt.show()

        
fig, axs = plt.subplots(1,1)
axs.plot(t[0:numsampleToPlot],my_yT[0:numsampleToPlot],'b+')
axs.plot(t10[0:numsampleToPlot*10],my_yT10[0:numsampleToPlot*10],'r')
axs.grid()
axs.set_xlabel('time(sec)')
axs.set_ylabel('y(t)')
plt.show()


"""
Lets  generate the Fourier Analysis Coefficients C_k
 and plot it as combine cosine Fourier Series representations
"""

YT  = np.fft.fft(my_yT)
magYT = abs(YT)/len(YT)


""" Lets generate the cosine combine representation for the FS of yt
    we need to take the  coeff 0 , and 2*coeff(1:N/2)
"""

HalfN = int(np.floor(len(YT)/2))
magYTRHS    = 2.*magYT[0:HalfN]; 
magYTRHS[0] = magYT[0];  # the first element is DC -> should not be doubled 
phaseYT = fn_PostProcessPhase(YT);
fig, axs = plt.subplots(2,1)
axs[0].stem(magYTRHS[0:highestFreq])
axs[1].grid()
axs[0].set_title('Combined Trigo Form FS')
axs[0].set_xlabel('frequency')
axs[0].set_ylabel('mag')
axs[1].stem(phaseYT[0:highestFreq])
axs[1].grid()
axs[1].set_xlabel('frequency')
axs[1].set_ylabel('Phase(rad)')

for idx,myOneFSCoeff in enumerate(tstSigFS):
     myA = myOneFSCoeff[0]
     myFreq = myOneFSCoeff[1]
     myPhase = myOneFSCoeff[2]
     
     axs[0].plot(myFreq, myA,'r+')
     axs[1].plot(myFreq, myPhase,'r+')
         
plt.show()




""" Lets plot the FS of yt as complex exponential double sided.
"""

magYTcentered = np.fft.fftshift(magYT);
phaseYTcentered = np.fft.fftshift(phaseYT);
freq = np.fft.fftfreq(len(magYT), d=1/Fs)
fshift= np.fft.fftshift(freq)
fshiftStartIdx = int((np.ceil(N/2))-highestFreq+2);
fshiftEndIdx = int((np.ceil(N/2))+highestFreq);


fig, axs = plt.subplots(2,1)
axs[0].stem(fshift[fshiftStartIdx:fshiftEndIdx],magYTcentered[fshiftStartIdx:fshiftEndIdx]); 
axs[0].grid(); 
axs[0].set_title('Complex Exponential Form FS')
axs[0].set_xlabel('frequency')
axs[0].set_ylabel('mag')

axs[1].stem(fshift[fshiftStartIdx:fshiftEndIdx],phaseYTcentered[fshiftStartIdx:fshiftEndIdx]); 
axs[1].grid()
axs[1].set_xlabel('frequency')
axs[1].set_ylabel('Phase(rad)')
plt.show()


