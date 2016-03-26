%this is script to test the correctness of the sine swept function
clc; clear all; close all;

initial_phase = 0;
starting_freq = 1;
final_freq = 200;
epsilon = 0.25;
amplitude = 1;

syms t;
%linear chirp, rate of frequency increase or chirp rate
k = (final_freq - starting_freq)/epsilon;
input_chirp = amplitude*sin(initial_phase + 2*pi*(starting_freq*t + (k*t^2)/2));

t = 0:1/(4*pi*final_freq):epsilon;
plot(t,subs(input_chirp));