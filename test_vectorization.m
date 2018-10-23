close all;
clear all;

t = 0:1/100:2*pi;
y = sin(t);
figure(1);
plot(t,y); hold on;
foundIdx = abs(y)>0.8;
% purpose - quickly find values |y| > 0.8 and set them to zero
y2 = y;
y2(abs(y)>0.8) = 1
plot(t,y2,'g+');
