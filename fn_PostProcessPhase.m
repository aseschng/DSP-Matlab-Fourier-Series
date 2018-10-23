function processedPhaseYT = fn_PostProcessPhase(YT)

    magYT = abs(YT)/length(YT);
    maxYT = max(magYT);
    processedPhaseYT = angle(YT).*(magYT>(1/1000)*maxYT);
    % making the phase of the +ve frequency near pi to become -pi
    % making the phase of the -ve frequency near pi to become +pi

    halfN = ceil(length(YT)/2);
    A = processedPhaseYT;
    Arange =[1:length(YT)];

    rangeLeftA = (abs((abs(A)-pi))<1e-6) & (Arange<halfN);
    A(rangeLeftA) = -pi;
    rangeRightA = (abs(abs(A)-pi)<1e-6) & (Arange>=halfN);
    A(rangeRightA) = +pi;
    
    rangeA = (abs(A-0)<1e-6);
    A(rangeA) = 0;
    processedPhaseYT = A;
    
end

