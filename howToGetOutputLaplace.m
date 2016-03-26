clc; clear all; close all;
%testing the way to get output
syms t s;
%first to specify the inputs
Total_Impulse = 1;
epsilon = 0.5;          %sort of defines the hardness of the hammer
input_square = (Total_Impulse/epsilon)*(heaviside(t) - heaviside(t - epsilon));

Amplitude_force = 1;
epsilon = 0.5;
input_sine_pulse = (Amplitude_force*sin(pi*t/epsilon))*(heaviside(t) - heaviside(t - epsilon));

initial_phase = 0;
starting_freq = 1;
final_freq = 20;
epsilon = 2;
amplitude = 1;
k = (final_freq - starting_freq)/epsilon;
input_chirp = amplitude*sin(initial_phase + 2*pi*(starting_freq*t + (k*t^2)/2));

%single order spring, mass, damper system
natural_frequency = 10;                %100 hertz
damped_natural_frequency = 9;          %90  hertz
damping_ratio = sqrt(1 - damped_natural_frequency^2/natural_frequency^2);
wd = 2*pi*damped_natural_frequency;
k = 1; %1 kn/metre, spring stiffness;
wn = 2*pi*natural_frequency;
mass = k/(wn^2);
c = 2*damping_ratio*sqrt(k*mass);
FRF = 1/(mass*s^2 + c*s + k);

%transform input to laplace domain
Input_square = laplace(input_square, t, s);
Input_sine_pulse = laplace(input_sine_pulse,t,s);
Input_chirp = laplace(input_chirp, t, s);

%get outputs in laplace domain then in time domain
Output_square = Input_square*FRF;
Ouput_sine_pulse = Input_sine_pulse*FRF;
Output_chirp = Input_chirp*FRF;

output_square = ilaplace(Output_square, s, t);
output_sine_pulse = ilaplace(Ouput_sine_pulse, s, t);
output_chirp = ilaplace(Output_chirp, s, t);


%specify parameters in frequency domain to perform discretization
Num_samples = 512;
Max_sampling_frequency = 2*pi*30; %3 times the frequency of the input
Duration = Num_samples/Max_sampling_frequency;
t = linspace(0, Duration, Num_samples);

%plot stuff to verify
figure(1);
plot(t, subs(output_square),'b'); hold on;
plot(t, subs(input_square), 'r');
legend('output', 'input');

figure(2);
plot(t, subs(output_sine_pulse),'b'); hold on;
plot(t, subs(input_sine_pulse), 'r');
legend('output', 'input');

% figure(3);
% plot(t, subs(output_chirp),'b'); hold on;
% plot(t, subs(input_chirp), 'r');
% legend('output', 'input');

figure(3);
syms w f;
s = -w*1i;
% w = 2*pi*f;
FRF = subs(FRF,s);
% FRF = subs(FRF,w);
w = linspace(0, Max_sampling_frequency, Num_samples);
plot(w, subs(abs(FRF)));




