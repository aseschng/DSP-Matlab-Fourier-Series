% this function use close form solution to get the various signal's FS
% coefficients as combined trigo for of th Fourier series, i.e sum of
% scaled cosine at various harmonics and its phase.
% you should/can modify the code for DIY part only as the others are
% correct close form for Square, Triangle and Saw.

function [myA,myF,myPhi,K] = fn_getVariousSignals_FS_Coeff(typeSignal, numK, N, myFundamentalFreq, optionCell)
% The variables are:
% para 1= 'DIY', 'Square', 'Saw', 'Triangle'
% Para 2= number of sinusoids to use (includes DC) 
% Para 3= number of samples in 1 sec (sampling freq)
% Para 4= fundamental freqeuncy of the periodic signal to construct
% Para 5= optionCell for various plots, e.g, DutyCycle for 'Square',
%          'Ascending' == {1,0} for Saw


w0  = 2*pi/N;
% lets specify the sinusoids with respect to w0.
% The equation below is from Math2, slide pg 64, Fourier Series Slides
% y(t) = 10 +3cos(w_0 t) + 5cos(2.5w_0 t + phi/6) + 4sin(3w_0 t)
if (strcmp(typeSignal,'DIY') == 1)
    K = 4;  % Lets have three sinusoid
    myA(1) = 10;   myF(1)   = 0*myFundamentalFreq;  myPhi(1) =  0;   
    myA(2) = 3;   myF(2)   = myFundamentalFreq;    myPhi(2) =  0;   
    myA(3) = 5;   myF(3)   = 2*myFundamentalFreq;  myPhi(3)   = +pi/6;  
    myA(4) = 4;   myF(4)   = 3*myFundamentalFreq;  myPhi(4)   = -pi/2;  
end  % of  DIY


if (strcmp(typeSignal,'Square') == 1)
% defining the coefficients of a square wave
% a[n] = (2*A/(n*pi))*sin(n*pi/2)
% where n = 1,2,3,4...
    K = numK;  
    A = 1;
    tmpCell     = optionCell{1}(1);  % This cell should be 'DutyCycle'
    tmpValCell  = optionCell{1}(2);  % This cell should be a value between 0->1.0
    DutyCycle   = tmpValCell{1}(1);
    
    % The dc term
    myA(1)  = A*DutyCycle;  myF(1) = 0;   myPhi(1) = 0;
    % Must re-calculate the DC coefficient with duty cycle.
    
    for (k=1:K-1)
        myA(k+1)   = (2*A/(k*pi))*sin(k*pi*DutyCycle);
        myF(k+1)   = k*myFundamentalFreq;  
        myPhi(k+1)  = 0;  

        if ( myA(k+1) < 0)
            myA(k+1) = abs(myA(k+1));
            myPhi(k+1) = -pi;
            if (abs(myA(k+1))<1e-6)
                myPhi(k+1) = 0;
            end  % of abs 
        end  % of myA
    end  % of for k
end  % of Square


if (strcmp(typeSignal,'Triangle') == 1)
% defining the coefficients of a saw Triagle
% a[n] = -(8*A)/(pi^2)*n^2)*((-1)^(n-1)/2)
% where n = 1,3,5  (odd values only)...

    K = numK;  A = 1;
    myA(1)  = A;  myF(1) = 0;   myPhi(1) = 0;
    for (k=1:K-1)
        if (mod(k,2) == 1)
            myA(k+1)   = ((8*A)/(pi*pi*k*k));
        else
            myA(k+1)   = 0;
    end
        myF(k+1)   = k*myFundamentalFreq;  
        myPhi(k+1) =  0;  
    end  % end of for k
end  % of Triangle


if (strcmp(typeSignal,'Saw') == 1)
% defining the coefficients of a saw Toothwave
    K = numK;  A = 1;
    myA(1)  = A/2;  myF(1) = 0;   myPhi(1) = 0;
    tmpCell = optionCell{1}(1);  % This cell should be 'Ascending'
    tmpValCell  = optionCell{1}(2);  % This cell should be a value between 0-> descending, 1-> ascending
    AscendingFlag  = (tmpValCell{1}(1) == 1);

    for (k=1:K-1)
        if (AscendingFlag ==1)
            myA(k+1)   = (2*A)/(k*2*pi);
            myPhi(k+1) =  +pi/2;  
        else
            myA(k+1)   = (2*A)/(k*2*pi);
            myPhi(k+1) =  -pi/2;  
        end
        myF(k+1)   = k*myFundamentalFreq;  

    end  % of for k
end  % of SAW

end  % of fn_getVariousSignals_FS

