%fourier transform trial script to check basic concepts that have been read
clc; clear all; close all;
Amp = 1;
freq = 12000;  %frequency of the sine wave
N = 256;
TimeLength = 1; %2 seconds worth of signal to be sampled
fsHz = N/TimeLength;
dt = TimeLength/N;

% sine = Amp*sin(2*pi*freq*(0:dt:TimeLength));
t = 0:dt:TimeLength;
% sine = Amp*(1 -cos(2*pi*5*t./3)).*cos(2*pi*5*t).*(t >= 0 & t <= 3/5);
sine = cos(40*pi*t) + 5*sin(pi/3 + 10*pi*t);

transform = fft(sine,N)/(N);
magTransform = abs(transform);
max(magTransform)
find(magTransform == max(magTransform))
faxis = linspace(0,fsHz/1000,N);
plot(faxis,(magTransform));
% axis([-40 40 0 0.6])
xlabel('Frequency (KHz)')