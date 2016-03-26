%this script is built to test if the inverse obtained by using the real
%observation matrix matches up with the initially sampled points

%define the number of sampled points of the function
Num_sample_points = 6;
%define the duration of the time
t_start = 0;
t_end = 0.5;

%discrete time sample points
k = linspace(t_start, t_end, Num_sample_points);

%define the function that is being used for this test
 f = 5 + 2*cos(pi*k/2 - pi/2) + 3*cos(pi*k);
 
 %get the fourier transform of the signal
 Y = dftmtx(Num_sample_points)*transpose(f);
 
 %separate the fourier transform into real and imaginary part in matrix
 counter = 1;
 X = zeros(Num_sample_points,1);
 for i=1:Num_sample_points
     X(counter) = real(Y(i));
     counter  = counter +1;
     X(counter) = imag(Y(i));
     counter = counter+1;
 end
     
 %get the observation matrix
 [H] = GetObservationMatrix(Num_sample_points)*Num_sample_points;
 
 %get the inverse transform using the observation matrix
 A = H*X;
 
 %check if A is the same as f
 tolerance = 0.000001;
 if(abs(sum(A - transpose(f))) < tolerance)
     fprintf('Everything works fine\n');
 else
     fprintf('Incorrect maths somewhere\n');
 end
 
 %test the full observation matrix
 sampled_number = 1;
 [H_complete] = CompleteObservationMatrix(H(sampled_number,:), Num_sample_points);
 
 %test the exponential co-variance
 ini_std = 10;
 dimension = Num_sample_points*2;
 psi = 0;
 scale = ini_std*ini_std;
 phi = 10;
 [Covariance] = ExponentialCoVariance1(ini_std,dimension, psi, scale, phi);
 
 
 