%This is the attempt to rectify the concept in KalmanFrequencyTransform by
%dividing the random variables to express the real and imaginary part of
%the frequency 

%The concept is probably still wrong. This should eventually lead on to 
%complex Kalman Filter (widely linear equations).
function KalmanFrequencyTransform_Intermediate()
%assumption is that the input signal is band-limited signal and can be
%adequately and completely represented by k bins in frequency domain
%number of samples = 2*B*L; fs (sampling frequency) = 2B; 
%fourier transform of signal Ff(s) is zero (negligible energy, beyond ability
%to measure) outside 0<s<2B range

%define variables for sampling the signal
Num_sample_points = 32;               %power of two here for easy FFT with no padding
time_start = 0;                        %start from 0 second
time_length = 4;                       %the end time of the signal for sampling (seconds), (L)
t = linspace(time_start, time_start + time_length, Num_sample_points+1);
t = t(1:end-1);
sampling_freq = t(2) - t(1);

% define input signal to be tested (any)
amplitude = 1;          %amplitude of the pulse
F = 5;                  %dominant frequency in the pulse
N = 3;                 %number of cycles (amount of "ringing") in the pulse and hence the bandwidth
input_signal = pulse_ref(amplitude,F,N,t);  %already sampled with t (discrete signal)

% max_height = 1;           %maximum height of the triangle
% t_mid = 0.75;             %time between +slope and -slope  
% t_end = 2;                %time after which y is zero  
% input_signal = triangle_ref(max_height, t_mid, t_end, t); %already sampled with t (discrete signal)

%define initial belief about our random variables, mean assumed to be zero (no contribution in frequency components)
%here the random variable is 2nx1 to encompass real and imaginary part of the complex variable
Xe = zeros(2*Num_sample_points, 1);   %number of frequency bins = number of sample points
ini_std = 50;
psi = 0;
scale = ini_std*ini_std;
phi = 10;
P = ExponentialCoVariance(ini_std,2*Num_sample_points, psi, scale, phi);  %dependence of bins that are relatively close
A = ones(1,Num_sample_points*2)*ini_std*ini_std;
P = diag(A);
%the noise associated with the measurement (polluting the time domain) is
%assumed to be zero mean white gaussian noise with standard deviation given below
std_sampled_signal = 0.3;   %in terms of amplitude

num_observations = Num_sample_points*200;  %for updating the initial belief
observations_per_sample = zeros(1,Num_sample_points);

%initialize plots of all figures
[plot_handles] = Iniplots(t,input_signal, sampling_freq, Num_sample_points, Xe, std_sampled_signal);

f = linspace(0, 1/sampling_freq, Num_sample_points+1);
f = f(1:end-1);

%initialize some variables for EKF update
% noisy_observation = zeros(Num_sample_points+1,1);
% v = ones(1,Num_sample_points+1).*std_sampled_signal.*std_sampled_signal.*4;
% R = diag(v);
% A = zeros(1,Num_sample_points/2-2);
% IDFT = zeros(Num_sample_points,Num_sample_points*2);

%don't use update for constraints
% noisy_observation = zeros(1,1);
R = std_sampled_signal*std_sampled_signal;

%define matrix for observation (real) as the inverse discrete fourier transform matrix
%keeping the constraint that the final signal is real in mind
IDFT= GetObservationMatrix(Num_sample_points);

%performing update steps Kalman filter
for i=1:num_observations
    %pick a random observation from the
    sampled_number = randi(Num_sample_points);
    
    %pollute sampled input signal at the sample number
    noisy_observation = Pollute(input_signal(sampled_number), std_sampled_signal);
    
    %get the observation matrix for this instance of update
%     H = CompleteObservationMatrix(IDFT(sampled_number,:), Num_sample_points);
    H = IDFT(sampled_number,:);

    %update step of kalman filter
    expected_measurement = H*Xe;
    %innovation, residual matrix
    z = noisy_observation - expected_measurement;
  
    S = R + H*P*H';
    iS = inv(S);
    
    K = P*H'*iS; %kalman gain
    
    Xe = Xe + K*z;                  %updated expected value
    P = P - P*H'*iS*H*P;            %updated co-variance matrix
    
    
    %update the plots
%     set(plot_handles(2), 'xdata', t, 'ydata', input_signal);
%     set(plot_handles(5), 'xdata', t(sampled_number), 'ydata', input_signal(sampled_number), 'MarkerEdgeColor', 'g');
    set(plot_handles(5), 'xdata', t(sampled_number), 'ydata', noisy_observation);
    %reshape matrix to separate matrix out
    Xe_separate = reshape(Xe, 2, Num_sample_points);
    set(plot_handles(3), 'xdata', f , 'ydata', Xe_separate(1,:));
    
    observations_per_sample(sampled_number) = observations_per_sample(sampled_number) + 1;
    set(plot_handles(4), 'xdata',1:Num_sample_points , 'ydata', observations_per_sample);
    
    
    %put a pause to be able to see how the update is happening
    pause(0.01);
end
display(isreal(P));





function [plot_handles] = Iniplots(t, input_signal, sampling_freq, Num_sample_points, Xe, std_noise)
%figure 1 just plots the actual frequency components (fourier transform) of the input signal
yf1 = FourierT(input_signal,sampling_freq, Num_sample_points);
f = linspace(0, 1/sampling_freq, Num_sample_points + 1);
f = f(2:end);
figure(1);
plot_handles(1) = plot(f,(yf1));
hold on;
title('Fourier Transform of Input Signal without noise');
xlabel('Frequency bins (k)');
ylabel('Amplitude');

%figure 2 shows the point in the time domain that is being used currently for the observation
figure(2);
plot_handles(2) = plot(t, input_signal, 'ro');
hold on;
plot_handles(5) = plot(0,0, 'g.', 'MarkerSize', 15); 
set(plot_handles(2), 'MarkerFaceColor', 'r');
plot(t, input_signal, 'k');
title('Actual Input Signal (Impulse)');
legend('Sampled discrete points','Current point being used in update step', 'Input Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');

%figure 3 shows the current mean of all the frequency bins random variables after each update
figure(3); %here, you use the initialized mean (could be zero)
plot_handles(3) = plot(Xe);
hold on;
title('Estimated amplitude (by Kalman filter) at various frequency bins');
xlabel('Frequency bins (k)');
ylabel('Amplitude');

%figure 4 shows the graph depicting the number of times each particular observation has been picked for update step
figure(4);
plot_handles(4) = bar(0,0);
hold on;
title('Numer of times each Sampled point is used for Kalman update step');
xlabel('Index of Sampled signal');
ylabel('Number of times picked');

%figure 5 shows the frequency components if the input signal were polluted by noise
figure(5);
noisy_input = Pollute(input_signal, std_noise);
yf1 = FourierT(noisy_input,sampling_freq, Num_sample_points);
plot(f,(yf1));
hold on;
title('Fourier Transform of Input Signal with noise');
xlabel('Frequency bins (k)');
ylabel('Amplitude');

%adds noise to the sampled input signal
function [noisy_input_signal] = Pollute(input_signal, std_noise)
noiseInMeasurements = std_noise*randn(size(input_signal));
noisy_input_signal = input_signal + noiseInMeasurements;


%generate pulse input signal for the experiment
function y = pulse_ref(A,F,N,t)
y = A*(1 -cos(2*pi*F*t./N)).*cos(2*pi*F*t).*(t >= 0 & t <= N/F);

%generate triangular shaped signal for the experiment
function y = triangle_ref(max_height, t_mid, t_end, t)
%calculate some intermediate parameters
slope_down = -max_height/(t_end-t_mid);
slope_up = max_height/t_mid;
c = -slope_down *t_end;

y = slope_up*t.*(t<t_mid) + (slope_down*t + c).*(t>=t_mid & t<=t_end);

%perform fast fourier transform of signal with sampling period dt
function y = FourierT(x,dt,Num_sample_points)
y = fft(x,Num_sample_points);


%perform fast inverse fourier transform of signal with sampling period dt
function y = IFourierT(x,dt,Num_sample_points)
y = ifft(x,Num_sample_points)/dt;

