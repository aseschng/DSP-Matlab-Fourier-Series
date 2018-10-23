function [myPhi_delayed] = fnDelaySignal_FS(N, myF, myPhi, vDelayTimeNormalized)

numSampleDelayed = vDelayTimeNormalized*N;
myPhi_delayed = myPhi;
    for (i=1:length(myF))
        whichFreq = myF(i);
        myPhi_delayed(i) = myPhi(i) - ((2*pi/N)*whichFreq*numSampleDelayed);
    end
end

