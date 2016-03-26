clc; clear all; close all;
%testing the way to get output
syms t w;
%first to specify the input
% input = sin(2*pi*10*t); %+ cos(2*pi*15*t);
Total_Impulse = 1;
epsilon = 0.5;          %sort of defines the hardness of the hammer
input = (Total_Impulse/epsilon)*(heaviside(t) - heaviside(t - epsilon));

amplitude_force = 1;
input_sine = (amplitude_force*sin(pi*t/epsilon))*(heaviside(t) - heaviside(t - epsilon));

%specify frequency response you want to test, low pass filter
cut_off_freq = 10;
FRF = heaviside(w + cut_off_freq) - heaviside(w - cut_off_freq);

%single order spring, mass, damper system
natural_frequency = 10; %100 hertz
damped_natural_frequency = 9;
damping_ratio = sqrt(1 - damped_natural_frequency^2/natural_frequency^2);
wd = 2*pi*damped_natural_frequency;
k = 1; %10 n/metre, spring stiffness;
syms w fn ;
% w = 2*pi*f;
wn = 2*pi*fn;
FRF = 1/(k*(1 + 2*1i*damping_ratio*w/wn - w^2/wn^2));
fn = natural_frequency;
FRF = subs(FRF);

%specify parameters in frequency domain to perform discretization
Num_samples = 128;
Max_sampling_frequency = 30;
Duration = Num_samples/Max_sampling_frequency; %3 times the frequency of the input
t = linspace(0, Duration, Num_samples);
w = linspace(0, Max_sampling_frequency, Num_samples);
input_subs = subs(input,t);
input_sine_subs = subs(input_sine,t);
figure(1);
plot(t, input_subs);hold on;
plot(t, input_sine_subs);

%transform input to frequency domain
Input = dftmtx(Num_samples)*input_subs';
figure(2);
plot(w, abs(Input)); hold on;
Input_sine = dftmtx(Num_samples)*input_sine_subs';


FRF_subs = subs(FRF,w);
figure(5);
plot(w, abs(FRF_subs));

Output = Input.*FRF_subs';
Output_sine = Input_sine.*FRF_subs';
figure(3);
plot(w, abs(Output));

figure(6);
% plot(w, abs(Output_sine)); hold on;
plot(w, abs(Input_sine),'r'); 

figure(4);
output = conj(dftmtx(Num_samples)/Num_samples)*Output;
output_sine = conj(dftmtx(Num_samples)/Num_samples)*Output_sine;
plot(t, real(output));hold on;
plot(t, real(output_sine));
