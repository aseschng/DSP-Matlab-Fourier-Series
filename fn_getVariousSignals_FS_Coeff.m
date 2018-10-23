function [myA,myF,myPhi,K] = fn_getVariousSignals_FS_Coeff(typeSignal, numK, N, myFundamentalFreq, optionCell)
% The variables are:
% para 1= 'DIY', 'Square', 'Saw'
% Para 2= number of sinusoids to use (includes DC) 
% Para 3= number of samples in 1 sec (sampling freq)
% Para 4= fundamental freqeuncy of the periodic signal to construct
% Para 5= optionCell for various plots, e.g, DutyCycle for 'Square',
%          'Ascending' == {1,0} for Saw


w0  = 2*pi/N;
% lets specify the sinusoids with respect to w0.
if (strcmp(typeSignal,'DIY') == 1)
    K = 3;  % Lets have three sinusoid
    myA(1) = 3;   myF(1)   = 0*myFundamentalFreq;  myPhi(1) =  0;   
    myA(2) = 2;   myF(2)   = myFundamentalFreq;    myPhi(2) =  +pi/3;   
    myA(3) = 1;   myF(3)   = 2.5*myFundamentalFreq;  myPhi(3)   = -pi/2;  
end


if (strcmp(typeSignal,'Square') == 1)
% defining the coefficients of a square wave
% a[n] = (2*A/(n*pi))*sin(n*pi/2)
% where n = 1,2,3,4...
    K = numK;  
    A = 1;
    tmpCell = optionCell{1}(1);  % This cell should be 'DutyCycle'
    tmpValCell  = optionCell{1}(2);  % This cell should be a value between 0->1.0
    DutyCycle = tmpValCell{1}(1);
    
    % The dc term
    myA(1)  = A/2;  myF(1) = 0;   myPhi(1) = 0;
    for (k=1:K-1)
        myA(k+1)   = (2*A/(k*pi))*sin(k*pi*DutyCycle);
        myF(k+1)   = k*myFundamentalFreq;  
        myPhi(k+1)  = 0;  

        if ( myA(k+1) < 0)
            myA(k+1) = abs(myA(k+1));
            myPhi(k+1) = -pi;
            if (abs(myA(k+1))<1e-6)
                myPhi(k+1) = 0;
            end
            
        end
    end
end


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
    end
end


if (strcmp(typeSignal,'Saw') == 1)
% defining the coefficients of a saw Toothwave
% a[n] = -(2*A/(n*pi))*(-1)^n
% where n = 1,2,3,4...
    K = numK;  A = 1;
    myA(1)  = 3;  myF(1) = 0;   myPhi(1) = 0;
    tmpCell = optionCell{1}(1);  % This cell should be 'Ascending'
    tmpValCell  = optionCell{1}(2);  % This cell should be a value between 0-> descending, 1-> ascending
    AscendingFlag  = (tmpValCell{1}(1) == 1);

    for (k=1:K-1)
        if (AscendingFlag ==1)
            myA(k+1)   = (2*A)/(k*pi);
            myPhi(k+1) = +pi/2;  
        else
            myA(k+1)   = (2*A)/(k*pi);
            myPhi(k+1) =  -pi/2;  
        end
        myF(k+1)   = k*myFundamentalFreq;  
    end
end



end

