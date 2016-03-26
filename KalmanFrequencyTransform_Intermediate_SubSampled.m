%This is the attempt to rectify the concept in KalmanFrequencyTransform by
%dividing the random variables to express the real and imaginary part of
%the frequency 

%The concept is probably still wrong. This should eventually lead on to 
%complex Kalman Filter (widely linear equations).
function KalmanFrequencyTransform_Intermediate_SubSampled()
%assumption is that the input signal is band-limited signal and can be
%adequately and completely represented by k bins in frequency domain
%number of samples = 2*B*L; fs (sampling frequency) = 2B; 
%fourier transform of signal Ff(s) is zero (negligible energy, beyond ability
%to measure) outside 0<s<2B range

%define variables for sampling the signal
Num_sample_points = 128;               %power of two here for easy FFT with no padding
Max_duration = 1;
syms t;
Max_sampling_freq = Num_sample_points/Max_duration;
division_factors = [1];
% define input signal to be tested (any)
amplitude = 1;          %amplitude of the pulse
F = 5;                  %dominant frequency in the pulse
N = 2;                 %number of cycles (amount of "ringing") in the pulse and hence the bandwidth
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
std_sampled_signal = 0;   %in terms of amplitude

% num_observations = Num_sample_points*20;  %for updating the initial belief
% observations_per_sample = zeros(1,Num_sample_points);

%initialize plots of all figures
[plot_handles] = Iniplots(input_signal, Max_sampling_freq,Max_duration, Num_sample_points, Xe, std_sampled_signal);



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

num_signal = 10;


max_index = SubSampleMatrix(Max_duration, Max_duration/division_factors(1), Max_sampling_freq);
max_samples = Max_duration*Max_sampling_freq/division_factors(1);


%performing update steps Kalman filter
for i=1:num_signal
      
    %sub_sample the signal
    [duration sampling_frequency] = SubSampler(Max_sampling_freq, Max_duration, division_factors, 1);
%     sampling_frequency = sampling_frequency/2;
    duration = duration/2;
%     sampling_frequency = Max_sampling_freq;
    num_samples = duration*sampling_frequency;
    t = linspace(0, duration, num_samples+1);
    t = t(2:end);
    Input_signal = subs(input_signal, t);
        
    %pollute sampled input signal at the sample number
    noisy_observation = Pollute(Input_signal, std_sampled_signal);
    
%     %get observation matrix according to the new number of sample points
%     IDFT= GetObservationMatrix(num_samples);

    %adapt IDFT for the new number of samples
    IDFT_adapted = IDFT.*Num_sample_points/num_samples;
    %update the belieft of the system with each sampled point
    [index] = SubSampleMatrix(Max_duration, duration, sampling_frequency);
    for j=1:num_samples
        %get the observation matrix for this instance of update
        H = IDFT_adapted(j,index);
        
        %update step of kalman filter
        expected_measurement = H*Xe(index);
        %innovation, residual matrix
        z = noisy_observation(j) - expected_measurement;
        
        S = R + H*P(index,index)*H';
        iS = inv(S);
        
        K = P(index,index)*H'*iS; %kalman gain
        
        Xe(index) = Xe(index) + K*z;                  %updated expected value
        P(index,index) = P(index,index) - P(index,index)*H'*iS*H*P(index,index);            %updated co-variance matrix
    end

    
    f = linspace(0, sampling_frequency, num_samples + 1);
    f = f(2:end);
    
    %update the plots
    set(plot_handles(5), 'xdata', t, 'ydata', noisy_observation);
    %reshape matrix to separate matrix out
    Xe_separate = reshape(Xe(index), 2, num_samples);
    set(plot_handles(3), 'xdata', f , 'ydata', Xe_separate(1,:));
      
    
    %put a pause to be able to see how the update is happening
    pause(0.001);
end
display(isreal(P));

%sub sample the H-matrix accordingly
function [index_location] = SubSampleMatrix(Max_duration, Duration, Sampling_Frequency)
index_location = zeros(1, 2*Duration*Sampling_Frequency);
index_location(1:2) = [1 2];

index = 2;
count = 1;
factor = Max_duration/Duration;

while(index<2*Duration*Sampling_Frequency)
    index_location(index+1:index+2) = [count + 2*factor count+2*factor+1];
    count = count + 2*factor;
    index = index + 2;
end


%function sub-sampler to determine the appropriate duration and sampling frequency  
%of any particular input based on max sampling frequency, duration and
%allowable division 
function [duration sampling_frequency] = SubSampler(Max_sampling_frequency, Max_duration, Division_factors, num)
%generate sub-sampler parameters num times
duration = zeros(1,num);
sampling_frequency = zeros(1,num);
for i=1:num
    %choose 'm' as a random integer number between 1 and max_division factor
    m = randi([1 size(Division_factors,2)]);
    m = Division_factors(1,m);
    %In order for the sub-sampled signal to have fourier transform overlay on
    %top of the defined frequency bins 2 conditions must be fullfilled
    %First condition: new_duration = original_duration/m, 'm' being whole number
    duration(i) = Max_duration/m;
    
    %Second condition: new number of samples must be chosen such that N_new <= N_old/m
    max_new_num_samples = Max_sampling_frequency*Max_duration/m;
    max_new_num_samples = floor(max_new_num_samples);
    
    %randomize new number of samples to follow the second equatility by
    %specifying the lower and upper limit of random generation
    new_num_samples = randi([floor(Max_sampling_frequency*Max_duration/(Division_factors(end)+1)) max_new_num_samples]);
    
    new_num_samples = max_new_num_samples;
    
    sampling_frequency(i) = (new_num_samples/duration(i));
end


function [plot_handles] = Iniplots(input_signal, sampling_freq,max_duration, Num_sample_points, Xe, std_noise)
%figure 1 just plots the actual frequency components (fourier transform) of the input signal
% max_duration = max_duration/2;
% sampling_freq = sampling_freq/2;
% Num_sample_points = Num_sample_points/2;
t = linspace(0, max_duration, Num_sample_points+1);
t = t(2:end);
yf1 = FourierT(subs(input_signal,t),1/sampling_freq, Num_sample_points);
f = linspace(0, sampling_freq, Num_sample_points + 1);
f = f(2:end);
figure(1);clf;
plot_handles(1) = plot(f,real(yf1));
hold on;
title('Fourier Transform of Input Signal without noise');
xlabel('Frequency bins (k)');
ylabel('Amplitude');

%figure 2 shows the point in the time domain that is being used currently for the observation
figure(2);clf;
ezplot(input_signal, [0 max_duration]);
hold on;
plot_handles(5) = plot(0,0, 'g.', 'MarkerSize', 15); 
title('Actual Input Signal (Impulse)');
legend('Input Signal','Sampled discrete points');
xlabel('Time (seconds)');
ylabel('Amplitude');

%figure 3 shows the current mean of all the frequency bins random variables after each update
figure(3); clf;%here, you use the initialized mean (could be zero)
plot_handles(3) = plot(Xe);
hold on;
title('Estimated amplitude (by Kalman filter) at various frequency bins');
xlabel('Frequency bins (k)');
ylabel('Amplitude');

%figure 4 shows the graph depicting the number of times each particular observation has been picked for update step
% figure(4);
% plot_handles(4) = bar(0,0);
% hold on;
% title('Numer of times each Sampled point is used for Kalman update step');
% xlabel('Index of Sampled signal');
% ylabel('Number of times picked');

%figure 5 shows the frequency components if the input signal were polluted by noise
figure(5);clf;
noisy_input = Pollute(subs(input_signal,t), std_noise);
yf1 = FourierT(noisy_input,1/sampling_freq, Num_sample_points);
plot(f,real(yf1));
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
y = A*(1 -cos(2*pi*F*t./N)).*cos(2*pi*F*t).*(heaviside(t) - heaviside(t-N/F));

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

