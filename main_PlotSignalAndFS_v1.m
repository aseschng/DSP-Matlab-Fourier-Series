% Lets plot some periodic wavefiles in time,
% and then generate its Fourier representations.
% you can select by changing the variable option = ???;

% Author: Chng Eng Siong
% Date: 19 Oct 2018

clear all
close all

% some basic initialization
colormap ='rgbckrgbckrgbck';

%%%%%%%%%%%%%%%%%%%
% Lets generate some periodic signals, 
%   we now have  'DIY', 'Saw (ascending/descencing)', 'Square (with different duty cycle)'
%   The following routines is to get the FS coefficients ONLY
%%%%%%%%%%%%%%%%%%%

option = 1;
switch option
    case 1
        strSignalToGenerate = 'DIY'; 

    case 2
        strSignalToGenerate = 'Square'; 
        
    case 3
        strSignalToGenerate = 'Triangle'

    case 4
        strSignalToGenerate = 'Saw'; 
        
   otherwise
        disp('Choose : (1) DIY, (2) Square, (3) Triangle,  (4) Saw')
end



% The sampling rate is Fs
Fs  = 8000;
N   = Fs;
% Lets generate 1 second worth of data, N = sample rate
% therefore, when we take FFT, the frequency bin will be 1Hz for each bin.
% in other words, the fundamental frequency (digital be) exactly 2\pi/N
% therefore the resolution at each bin is 1 Hz! since we generate N samples
% for 1 second analysis window.
w0  = 2*pi/N;
myFundamentalFreq = 1000;  % we will generate periodic signals with this fundamental freqeuncies
K   = 10;  % number of sinusoids to generate
vDelayTimeNormalized = 0.0*(1/myFundamentalFreq);  % how much of 1 period of fundamental



switch strSignalToGenerate
    
    case 'DIY'
        [myA,myF,myPhi,K] = fn_getVariousSignals_FS_Coeff('DIY',4,N,myFundamentalFreq);

    case 'Square'
        % option For Square an array of cell, string and value pair, {{'DutyCycle',0..1 range}}
        [myA,myF,myPhi,K] = fn_getVariousSignals_FS_Coeff('Square',K,N,myFundamentalFreq,{{'DutyCycle',0.5}});

    case 'Saw'
        % optionStr {{'Ascending',0 or 1}} -> to plot ascending or descending saw
        [myA,myF,myPhi,K] = fn_getVariousSignals_FS_Coeff('Saw',K,N,myFundamentalFreq,{{'Ascending',1}});

        
    case 'Triangle'
        [myA,myF,myPhi,K] = fn_getVariousSignals_FS_Coeff('Triangle',K,N,myFundamentalFreq);
        
        
    otherwise
        disp('Choose : DIY, Square or Saw only')

end    



%%%%%%%%%%
%%%%%%%%%
% Lets generate 1 second worth of data in Fs sampling freqeuncy
% using the coefficients
%%%%%%%%%%%%%%%%%%%
if(vDelayTimeNormalized ~= 0)
    myPhi_orig    = myPhi;
    myPhi_delayed = fnDelaySignal_FS(N, myF, myPhi, vDelayTimeNormalized);
    myPhi         = myPhi_delayed;
    [t, my_y_orig, my_yT_orig] = fn_genTimeSignalFrom_FSCoeff(myA, myF, myPhi_orig, K, Fs);
    [t_delayed, my_y_delayed, my_yT_delayed] = fn_genTimeSignalFrom_FSCoeff(myA, myF, myPhi_delayed, K, Fs);
    
    figure(5);  
    numCycle=3;
    if (myF(1) ~= 0)
        numSampleToPlot = ceil(numCycle*(1/myF(1))*Fs);
    else
        numSampleToPlot = ceil(numCycle*(1/myF(2))*Fs);
    end
    st = 1; se = numSampleToPlot;
    plot(t(st:se),my_yT_orig(st:se),'b'); hold on;
    plot(t(st:se),my_yT_delayed(st:se),'r');
    title('plotting the original sequence (blue) vs delayed seq (red)');
    xlabel('time (sec)'); grid on; hold off;
end


%%%%%%%%%%%%%%%%%%%
% Lets generate 1 second worth of data in 10*Fs sampling freqeuncy
% This is to generate the continuous visualzation of the signal
%%%%%%%%%%%%%%%%%%%
% redo AGAIN the my_yT so that we are using the actual delayed myPhi!!!
[t, my_y, my_yT]       = fn_genTimeSignalFrom_FSCoeff(myA, myF, myPhi, K, Fs);
[t10, my_y10, my_yT10] = fn_genTimeSignalFrom_FSCoeff(myA, myF, myPhi, K, Fs*10);


%%%%%%%%%%%%%%%%%%%
% Lets plot the time domain representation of each sinusoid
% thats going to add up
%%%%%%%%%%%%%%%%%%%
numCycle = 3;
if (myF(1) ~= 0)
    numSampleToPlot = ceil(numCycle*(1/myF(1))*Fs);
else
    numSampleToPlot = ceil(numCycle*(1/myF(2))*Fs);
end
st = 1; se = numSampleToPlot;
figure(1);
tt = sprintf(' Time domain representation of y : numSinusoid = %d',K);
title(tt);
xlabel('time (sec)');
ylabel('y_k(t)'); hold on;
for (k=1:K)
    if (myPhi(k) >= 0)
        tmpLegend = sprintf('y%d(t) = %3.2f cos(2pi %d + %3.2f)',k, myA(k), myF(k), myPhi(k));
    else
        tmpLegend = sprintf('y%d(t) = %3.2f cos(2pi %d  %3.2f)',k, myA(k), myF(k), myPhi(k));
    end
    
    tt=sprintf('%c%c',colormap(mod(k,length(colormap))+1),char('o'));
%    plot(t(st:se),my_y(k,st:se),tt); 
    tt=sprintf('%c%c',colormap(mod(k,length(colormap))+1),char('-'));
    plot(t10(st:se*10),my_y10(k,st:se*10),tt);
    opLegend{k} = tmpLegend; hold on; grid on;
end
legend(opLegend);


%%%%%%%%%%%%%%%%%%%
% Lets plot the time domain representation of the summed sinusoid
%%%%%%%%%%%%%%%%%%%
figure(2);
hold on;
plot(t10(st:se*10),my_yT10(st:se*10),'g');
tt=sprintf(' Time domain representation yT(t) = sum of %d sinusoids',K);
title(tt);
xlabel('time (sec)');
ylabel('yT(t)'); grid on;



%%%%%%%%%%%%%%%%%%%
% Lets plot the fourier representation as cosines
% this is called the combined trigonimetric FS
%%%%%%%%%%%%%%%%%%
% Figure (4) -> the visualization of the FS coefficeients as sum of cosine
% + phase shift, hence
% Lets generate the Frequency domain representation ONLY the RHS
% hence it is plotting the \sum Magnitude*cosine(frequency + Phase)

YT  = fft(my_yT);
magYT = abs(YT)/N; 

% We extract ONLY the right hand side (+ve freqeuncies) hence must display
% with scaling * 2 -> this IS NOT the normal fourier transform, 
HalfN = floor((N-1)/2);
magYTRHS(1:HalfN) = 2.*magYT(1:HalfN); 
magYTRHS(1) = magYT(1);  % the first element is DC -> should not be doubled 
% see the table of combined Trigonometric

%phaseY = angle(YT);
phaseYT = fn_PostProcessPhase(YT);
phaseYTnormalized  = phaseYT/pi;

% plotting only for +ve frequency , showing A*cosine(frequency+phase);
figure;
fval = 0:N-1; % zero-centered frequency range
subplot(2,1,1); stem(fval(1:HalfN),magYTRHS); grid on;
title(' Magnitude plot of yT - as sum of cosine+phase shift');
xlabel('Frequency (Hz)'); grid on;

subplot(2,1,2); stem(fval(1:HalfN),phaseYT(1:HalfN)); grid on;
title(' Phase plot of yT');
ylabel('radian');
xlabel('Frequency (Hz)'); grid on;



%%%%%%%%%%%%%%%%%%%
%Lets plot the Fourier coefficients as complex exponentials
%%%%%%%%%%%%%%%%%%%
% Figure (5) -> the visualization of the FS coefficeients as sum of complex
% exponential
% Lets plot the Typical left and right vizualization in complex domain
YT  = fft(my_yT);
magYT = abs(YT)/N; 
magYTcentered = fftshift(magYT);
phaseYT = fn_PostProcessPhase(YT);
phaseYTcentered = fftshift(phaseYT);
fshift = (-N/2:N/2-1)*(Fs/N); % zero-centered frequency range

figure(4);
subplot(2,1,1); stem(fshift,magYTcentered); hold on; grid on; 
title(' Magnitude plot of yT - as sum of complex exponential');
xlabel('Frequency (Hz)'); 


subplot(2,1,2); stem(fshift,phaseYTcentered); hold on;grid on; 
title(' Phase plot of yT');
ylabel('radian');
xlabel('Frequency (Hz)'); grid on;



% Compare the FFT of the signal vs the sinsusoids that compose it!
% plotting only for +ve frequency , 
figure;
fval = 0:N-1; % zero-centered frequency range
subplot(2,1,1); stem(fval(1:HalfN),magYTRHS,'b-o');  hold on;
title(' Comparing FFT magnitude and the equation of sinusoids');
xlabel('Frequency (Hz)'); grid on;

subplot(2,1,2); stem(fval(1:HalfN),phaseYT(1:HalfN),'b-o'); hold on;
title(' Comparing FFT Phase and the equation of sinusoids');
ylabel('radian');
xlabel('Frequency (Hz)'); grid on;

for (k=1:K)
    tmpLegend = sprintf('y%d(t) = %3.2f cos(2pi %d + %3.2f)',k, myA(k), myF(k), myPhi(k));
    subplot(2,1,1); stem(myF(k),myA(k),'r+'); grid on; 
    subplot(2,1,2); stem(myF(k), myPhi(k),'r+'); grid on;
end


