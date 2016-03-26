%this is a dummy test to obtain the frequency response of a system using
%KF outlined in KalmanFrequencyTransform_Intermediate.m

%Here, the FRF is specified and known inputs and simulated outputs (both
%polluted with noise) are used to obtain FRF of the system.

%This concept should pave way to a full blown simulation program whereby user can:
%*modify FRF
%*specify characteristics of noises
%*specify input to be used for the experiment
%*modify frequency of sampling of different inputs
%*specify location of truncation of inputs
%*modify R matrix

%Assumption : Input and output polluted by the same noise since using the
%same measurement device, the truncation and sampling frequency on input and output is the same

function FRF_Kalman_Test()
%specify the FRF going to be used as the ground truth
[FRF SystemParams] = defineSystem();


%define the varuiables for frequency domain
Num_sample_points = 128;   %also the number of bins in frequency domain
Max_duration = 1.5;          %this is the maximum duration in the time domain, it defines the frequency resolution of fourier transform
Max_sampling_frequency = Num_sample_points/Max_duration;    %defines the maximum frequency represented by the fourier transform
Max_division_factor = [1];    %largest m, see function sub-sampler to understand more

%define the inputs keeping in mind the values defined above for the frequency domain
[input_square input_sine_pulse output_square output_sine_pulse]...
    = defineInputsOutputs(SystemParams);


%define initial belief and variables for kalman filter
Xe = zeros(2*Num_sample_points, 1);   %number of frequency bins = number of sample points
ini_std = 5;
% psi = 0;
% scale = ini_std*ini_std;  
% phi = 10;
% P = ExponentialCoVariance(ini_std,2*Num_sample_points, psi, scale, phi);  %dependence of bins that are relatively close
A = ones(1,Num_sample_points*2)*ini_std*ini_std;
P = diag(A);
%the noise associated with the measurement (polluting the time domain) is
%assumed to be zero mean white gaussian noise with standard deviation given below
std_sampled_signal = 0;   %in terms of amplitude, same for all inputs since using the same device

%initialize plots of all figures
[plot_handles Durations Sampling_Frequencies] = Iniplots(input_square, input_sine_pulse,...
    output_square, output_sine_pulse, FRF, Max_sampling_frequency, Max_duration, Max_division_factor);

Num_experiments = 20;

R = std_sampled_signal*std_sampled_signal;

IDFT= GetObservationMatrix(Num_sample_points);

%generate frequecy bins
f = linspace(0, Max_sampling_frequency, Max_sampling_frequency*Max_duration +1);
f = f(2:end);

%performing update step for the kalman filter
%performing update steps Kalman filter
for i=1:Num_experiments
    %based on the current sampling frequency and durations obtain the
    %discrete inputs and outputs
    Num_samples = Sampling_Frequencies.*Durations;
    t_square = linspace(0, Durations(1), Num_samples(1) + 1);
    t_square = t_square(2:end);
    t_sine_pulse = linspace(0, Durations(2), Num_samples(2) + 1);
    t_sine_pulse = t_sine_pulse(2:end);
    Input_square = subs(input_square.Data, t_square);
    Input_sine_pulse = subs(input_sine_pulse.Data, t_sine_pulse);
    Output_square = subs(output_square, t_square);
    Output_sine_pulse = subs(output_sine_pulse, t_sine_pulse);
    
    
    %pollute sampled input and output signals
    noisy_input_square = Pollute(Input_square, 0);
    noisy_output_square = Pollute(Output_square, std_sampled_signal);
    noisy_input_sine_pulse = Pollute(Input_sine_pulse, 0);
    noisy_output_sine_pulse = Pollute(Output_sine_pulse, std_sampled_signal);
    
    %perform the discete fourier transform of the inputs
    Noisy_input_square =  FourierT(noisy_input_square,double(Num_samples(1)));
    Noisy_input_square = AdaptMatrix(Noisy_input_square);
    Noisy_input_sine_pulse = FourierT(noisy_input_sine_pulse,double(Num_samples(2)));
    Noisy_input_sine_pulse = AdaptMatrix(Noisy_input_sine_pulse);
    
%     IDFT_adapted= IDFT.*(Num_sample_points/Num_samples(1));
%     
%     %have to select the frequency bands to update depending on duration and sampling frequency
%     [index] = SubSampleMatrix(Max_duration, Durations(1), Sampling_Frequencies(1));
%     %first update the belief of Frequency response function for the square inputs and outputs
%     for j=1:Num_samples(1)
%         %get the observation matrix for this instance of update
%         H = IDFT_adapted(j,index);
%               
%         %observation matrix has to be modified according to the input
% %         display([size(H) size(Noisy_input_square)]);
%         H = H.*Noisy_input_square;
%         
%         %update step of kalman filter
%         expected_measurement = H*(Xe(index));
%         %innovation, residual matrix
%         z = noisy_output_square(j) - expected_measurement;
%         
%         S = R + H*P(index,index)*H';
%         iS = inv(S);
%         
%         K = P(index,index)*H'*iS; %kalman gain
%         
%         Xe(index) = Xe(index) + K*z;                  %updated expected value
%         P(index,index) = P(index,index) - P(index,index)*H'*iS*H*P(index,index);            %updated co-variance matrix
%     end
    
    IDFT_adapted= IDFT.*(Num_sample_points/Num_samples(2));       
    %have to select the frequency bands to update depending on duration and sampling frequency
    [index] = SubSampleMatrix(Max_duration, Durations(2), Sampling_Frequencies(2));
    
    %update the belief with the sine_pulse inputs and outputs
    for j=1:Num_samples(2)
        %get the observation matrix for this instance of update
        H = IDFT_adapted(j,index);
        
        
        %observation matrix has to be modified according to the input
        H = H.*Noisy_input_sine_pulse;
        
        %update step of kalman filter
        expected_measurement = H*(Xe(index));
        %innovation, residual matrix
        z = noisy_output_sine_pulse(j) - expected_measurement;
        
        S = R + H*P(index,index)*H';
        iS = inv(S);
        
        K = P(index,index)*H'*iS; %kalman gain
        
        Xe(index) = Xe(index) + K*z;                  %updated expected value
        P(index,index) = P(index,index) - P(index,index)*H'*iS*H*P(index,index);            %updated co-variance matrix
    end
    %update the plots
    %for figure 1,
    set(plot_handles(1), 'xdata', t_square, 'ydata', Input_square);
    set(plot_handles(2), 'xdata', t_square, 'ydata', noisy_input_square);
    set(plot_handles(3), 'xdata', [Durations(1) Durations(1)],'ydata', [0 input_square.Total_impulse/input_square.Epsilon]);
    %for figure 2,
    set(plot_handles(4), 'xdata', t_sine_pulse, 'ydata', Input_sine_pulse);
    set(plot_handles(5), 'xdata', t_sine_pulse, 'ydata', noisy_input_sine_pulse);
    set(plot_handles(6), 'xdata', [Durations(2) Durations(2)],'ydata', [0 input_sine_pulse.Amplitude_force]);
    %for figure 4,
    set(plot_handles(7), 'xdata', t_square, 'ydata', Output_square);
    set(plot_handles(8), 'xdata', t_square, 'ydata', noisy_output_square);
    set(plot_handles(9), 'xdata', [Durations(1) Durations(1)],'ydata',[0 input_square.Total_impulse/input_square.Epsilon]);
    %for figure 5,
    set(plot_handles(10), 'xdata', t_sine_pulse, 'ydata', Output_sine_pulse);
    set(plot_handles(11), 'xdata', t_sine_pulse, 'ydata', noisy_output_sine_pulse);
    set(plot_handles(12), 'xdata', [Durations(2) Durations(2)],'ydata', [0 input_sine_pulse.Amplitude_force]);
    %for figure 8,
    %reshape matrix to separate matrix out
    Xe_separate = reshape(Xe, 2, Num_sample_points);
    set(plot_handles(13), 'xdata', f , 'ydata', Xe_separate(1,:));
    %get new sub_sampler for the inputs
    [Durations Sampling_Frequencies] = SubSampler(Max_sampling_frequency, Max_duration, Max_division_factor, 3);
    
    %put a pause to be able to see how the update is happening
    pause(0.1);
    fprintf('Iteration Num %d: \n', i);
   
end


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

%generate square impulse input signal
function [input_square] = Square_Impulse(Total_Impulse, epsilon)
syms t;    %defining the input function in terms of t to vary the sampling frequency later on
input_square = (Total_Impulse/epsilon)*(heaviside(t) - heaviside(t - epsilon));


%generate sine impulse signal
function [input_sine_pulse] = Sine_Impulse(Amplitude_Force, epsilon)
syms t;
input_sine_pulse = (Amplitude_Force*sin(pi*t/epsilon))*(heaviside(t) - heaviside(t - epsilon));


%generate linear sine sweep input signal
function [input_chirp] = Chirp(initial_phase, starting_freq, final_freq, duration, amplitude)
syms t;
%linear chirp, rate of frequency increase or chirp rate
k = (final_freq - starting_freq)/duration;
input_chirp = amplitude*sin(initial_phase + 2*pi*(starting_freq*t + (k*t^2)/2));

function [FRF SystemParams] = defineSystem()
%here, the system is a single order degree of freedom under damped spring, mass, damper system
natural_frequency = 10; %hertz
damped_natural_frequency = 9; %hertz
damping_ratio = sqrt(1 - damped_natural_frequency^2/natural_frequency^2);
spring_stiffness = 1; %1Kn/metre
syms w f fn;
w = 2*pi*f;
wn = 2*pi*fn;
FRF = 1/(spring_stiffness*(1 + 2*1i*damping_ratio*w/wn - w^2/wn^2));
fn = natural_frequency;
FRF = subs(FRF);

%put it the parameters of the single degree of freedom
SystemParams.fn = fn;
SystemParams.wn = 2*pi*fn;
SystemParams.k = spring_stiffness;
SystemParams.wd = 2*pi*damped_natural_frequency;
SystemParams.dr = damping_ratio;
SystemParams.m = spring_stiffness/(SystemParams.wn^2);
SystemParams.c = 2*damping_ratio*sqrt(spring_stiffness*SystemParams.m );


function [input_square input_sine_pulse output_square output_sine_pulse]...
    = defineInputsOutputs(SystemParams)

%first input is an impulse square wave of certain total impulse
input_square.Total_impulse = 1;
input_square.Epsilon = 1;          %sort of defines the hardness of the hammer
input_square.Data = Square_Impulse(input_square.Total_impulse, input_square.Epsilon);

%second input is a sinusoidal impulse of certain amplitude force
input_sine_pulse.Amplitude_force = 1;
input_sine_pulse.Epsilon = 1;
input_sine_pulse.Data = Sine_Impulse(input_sine_pulse.Amplitude_force, input_sine_pulse.Epsilon);

%third input is a chrip signal
input_chirp.Initial_phase = 0;
input_chirp.Starting_freq = 10;
input_chirp.Final_freq = 125;
input_chirp.Epsilon = 5;
input_chirp.Amplitude = 1;
input_chirp.Data = Chirp(input_chirp.Initial_phase, input_chirp.Starting_freq,...
    input_chirp.Final_freq, input_chirp.Epsilon, input_chirp.Amplitude);

%define FRF in laplace domain
syms t s;
FRF = 1/(SystemParams.m*s^2 + SystemParams.c*s + SystemParams.k);
%transform inputs to laplace domain
Input_square = laplace(input_square.Data, t, s);
Input_sine_pulse = laplace(input_sine_pulse.Data,t,s);

%get outputs in laplace domain then in time domain
Output_square = Input_square*FRF;
Ouput_sine_pulse = Input_sine_pulse*FRF;

output_square = ilaplace(Output_square, s, t);
output_sine_pulse = ilaplace(Ouput_sine_pulse, s, t);


%function to initialize plot variables
function [plot_handles Durations Sampling_Frequencies] = Iniplots(input_square, input_sine_pulse,...
    output_square, output_sine_pulse,FRF, Max_sampling_frequency, Max_duration, Max_division_factor)
%generate time samples with max_duration and max frequecy in mind
t = linspace(0, Max_duration, Max_sampling_frequency*Max_duration + 1);
t = t(2:end);
[Durations Sampling_Frequencies] = SubSampler(Max_sampling_frequency, Max_duration, Max_division_factor, 3);
%figures 1 to 3 would be the current input signals being used
figure(1);  %input_square wave
%new time based on it's sampling frequency
t_new = linspace(0, Durations(1), Sampling_Frequencies(1)*Durations(1) + 1);
t_new = t_new(2:end);
plot_handles(1) = plot(t_new, subs(input_square.Data, t_new), '.y'); hold on;
plot(t, subs(input_square.Data), 'b');
%the mark points for pollution with noise
plot_handles(2) = plot(0,0,'.r');
%plot the line showing the truncation line
plot_handles(3) = plot([Durations(1) Durations(1)],[0 input_square.Total_impulse/input_square.Epsilon], 'g'); 
title('Current impulse square input of system');
legend('Ideal Sampled Points', 'Input Signal', 'Actual Sampled Points (polled with noise)', 'Truncation line');

figure(2); %input_sine_pulse wave
%new time based on it's sampling frequency
t_new = linspace(0, Durations(2), Sampling_Frequencies(2)*Durations(2) + 1);
t_new = t_new(2:end);
plot_handles(4) = plot(t_new, subs(input_sine_pulse.Data, t_new), '.y'); hold on;
plot(t, subs(input_sine_pulse.Data), 'b');
%the mark points for pollution with noise
plot_handles(5) = plot(0,0,'.r');
%plot the line showing the truncation line
plot_handles(6) = plot([Durations(2) Durations(2)],[0 input_sine_pulse.Amplitude_force], 'g');
title('Current impulse sine input of system');
legend('Ideal Sampled Points', 'Input Signal', 'Actual Sampled Points (polled with noise)', 'Truncation line');

% figure(3); %input_chirp wave
% %new time based on it's sampling frequency
% t_new = linspace(0, Durations(3), Sampling_Frequencies(3)*Durations(3) + 1);
% t_new = t_new(2:end);
% plot(t_new, subs(input_chirp.Data,t_new), '.k'); hold on;
% plot(t, subs(input_chirp.Data), 'b');
% %the mark points for pollution with noise
% plot_handles(5) = plot(0,0,'.r');
% %plot the line showing the truncation line
% plot_handles(6) = plot([Durations(3) Durations(3)],[0 input_chirp.Amplitude], 'g');
% title('Current chirp input of system');
% legend('Ideal Sampled Points', 'Input Signal', 'Actual Sampled Points (polled with noise)', 'Truncation line');

%figures 4 to 6 would be the current output signals being obtained with same truncations as the corresponding inputs
figure(4); %output of input_square wave
%time based on sampling frequency of input
t_new = linspace(0, Durations(1), Sampling_Frequencies(1)*Durations(1) + 1);
t_new = t_new(2:end);
plot_handles(7) = plot(t_new, subs(output_square, t_new), '.y'); hold on;
plot(t, subs(output_square), 'b');
%the mark points for pollution with noise
plot_handles(8) = plot(0,0,'.r');
%plot the line showing the truncation line
plot_handles(9) = plot([Durations(1) Durations(1)],[0 input_square.Total_impulse/input_square.Epsilon], 'g'); 
title('Current impulse square response of system');
legend('Ideal Sampled Points', 'Output Signal', 'Actual Sampled Points (polled with noise)', 'Truncation line');

figure(5); %output of the input_sine_pulse wave
%time based on sampling frequency of input
t_new = linspace(0, Durations(2), Sampling_Frequencies(2)*Durations(2) + 1);
t_new = t_new(2:end);
plot_handles(10) = plot(t_new, subs(output_sine_pulse, t_new), '.y'); hold on;
plot(t, subs(output_sine_pulse), 'b');
%the mark points for pollution with noise
plot_handles(11) = plot(0,0,'.r');
%plot the line showing the truncation line
plot_handles(12) = plot([Durations(2) Durations(2)],[0 input_sine_pulse.Amplitude_force], 'g');
title('Current impulse sine response of system');
legend('Ideal Sampled Points', 'Output Signal', 'Actual Sampled Points (polled with noise)', 'Truncation line');

%figure 7 to represent the ground truth of the FRF being examiined
figure(7);
%generate frequecy bins
f = linspace(0, Max_sampling_frequency, Max_sampling_frequency*Max_duration +1);
f = f(2:end);
%plot absolute values of FRF
plot(f, subs(real(FRF)));

%figure 8 to represent the estimate of FRF at the moment
figure(8);
plot_handles(13) = plot(f, zeros(1,Max_sampling_frequency*Max_duration)); 


%function sub-sampler to determine the appropriate duration and sampling frequency  
%of any particular input based on max sampling frequency, duration and
%allowable division 
function [duration sampling_frequency] = SubSampler(Max_sampling_frequency, Max_duration, Division_factors, num)
duration = ones(1,num).*Max_duration;
sampling_frequency = ones(1,num).*Max_sampling_frequency;
return;
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
    
    sampling_frequency(i) = (new_num_samples/duration(i));
end

%perform fast fourier transform of signal with sampling period dt
function y = FourierT(x,Num_sample_points)
y = fft(x,Num_sample_points);

%adds noise to the sampled input signal
function [noisy_input_signal] = Pollute(input_signal, std_noise)
noiseInMeasurements = std_noise*randn(size(input_signal));
noisy_input_signal = input_signal + noiseInMeasurements;

%separate out the real and imaginary parts and align them
function [adapted] = AdaptMatrix(original)
dimension = size(original);
if(dimension(1) == 1)
    new_dimension = [1 dimension(2)*2];
    max = dimension(2);
else
    new_dimension = [dimension(1)*2 1];
    max = dimension(1);
end
adapted = zeros(new_dimension);
counter = 1;
for i = 1:max
    adapted(counter) = real(original(i));
    adapted(counter+1) = imag(original(i));
    counter = counter +2;
end

