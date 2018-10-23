function [t, my_y, my_yT] = fn_genTimeSignalFrom_FSCoeff(myA, myF, myPhi, K, Fs)

N  = Fs;   % we generate 1 second worth of data
w0 = 2*pi/N;
n=0:N-1;  % 1 seconds has 8000 samples;
t   = n.*1/Fs;
for (k=1:K)
  my_y(k,:) = myA(k)*cos(w0*myF(k)*n+myPhi(k));
end
my_yT  = sum(my_y);
end

