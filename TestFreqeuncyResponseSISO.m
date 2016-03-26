clc; clear all; close all;
%this is basically to test the frequency response of a single degree of
%freedom system with defined natural frequency, damping and stiffness

%define the FRF of SOF system 
natural_frequency = 10; %100 hertz
damped_natural_frequency = 9;
damping_ratio = sqrt(1 - damped_natural_frequency^2/natural_frequency^2);
wd = 2*pi*damped_natural_frequency;
k = 1; %1Kn/metre, spring stiffness;
syms w f fn s;
w = 2*pi*f;
wn = 2*pi*fn;
% s = w/wn;
FRF = 1/(k*(1 + 2*1i*damping_ratio*w/wn - w^2/wn^2));
fn = natural_frequency;
FRF = subs(FRF);

%define the FRF as a simple low pass filter
cut_off_freq = 50; %hertz
syms f;
% FRF = heaviside(f + cut_off_freq) - heaviside(f - cut_off_freq);

%define different inputs that you want to try on the system here defined in time domain
syms t;
Total_Impulse = 1;
epsilon = 0.5;          %sort of defines the hardness of the hammer
input_square = (Total_Impulse/epsilon)*(heaviside(t) - heaviside(t - epsilon));

Amplitude_Force = 1;
input_sine_pulse = (Amplitude_Force*sin(pi*t/epsilon))*(heaviside(t) - heaviside(t - epsilon));

initial_phase = 0;
starting_freq = 20;
final_freq = 100;
duration = 2;
amplitude = 1;
k = (final_freq - starting_freq)/duration;
input_chirp = amplitude*sin(initial_phase + 2*pi*(starting_freq*t + (k*t^2)/2));

%fourier transform the inputs to frequency domain
Input_square = fourier(input_square, t, f);
Input_sine_pulse = fourier(input_sine_pulse, t, f);
Input_chirp = fourier(input_chirp, t, f);

%obtain the output of the system in frequency domain
Output_square = Input_square*FRF;
Output_sine_pulse = Input_sine_pulse*FRF;
Output_chirp = Input_chirp*FRF;
Output_chirp = simplify(Output_chirp);

%obtain the outputs in time domain
Num_samples = 512;
Max_Duration = 1;
Maximum_sampling_frequency = Num_samples/Max_Duration;
f = linspace(0, Maximum_sampling_frequency, Num_samples +1);
f = f(2:end);
t = linspace(0, Max_Duration, Num_samples+1);
t = t(2:end);   
output_square = ifft(subs(Output_square,f), Num_samples)/Num_samples;

% output_square = ifourier(Output_square, f, t);
% output_sine_pulse = ifourier(Output_square, f, t);
% output_chirp = ifourier(Output_chirp, f, t);

%plot these things to see the response
plot(t, output_square);


