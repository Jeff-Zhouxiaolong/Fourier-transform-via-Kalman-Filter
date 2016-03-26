%testing the complete observation matrix in a small scale
function [H] = CompleteObservationMatrix(IDFT, Num_sample_points)
H(1,:) = IDFT;
%force the imaginary values to be zero here itself (probably not needed
%since update step does not involve solving for this parameter and there is
%no predictive step, so most likely stay as a prior zero anyway)
H(2,2) = 1;                             %DC component
H(3,Num_sample_points + 2) = 1;         %highest frequency component

counter = 1;
index = 4;
while (counter <Num_sample_points/2),
   % the real part of the n index and N-n index should be the same 
   H(index, counter*2+1) = 1;
   H(index, (counter+Num_sample_points/2)*2+1) = -1;
   %the imaginary part of the n index and N-n index should be conjugate
   index = index + 1;
   H(index, counter*2+2) = 1;
   H(index, (counter+Num_sample_points/2)*2+2) = 1;
   index = index + 1;
   
   %increment counter
   counter = counter + 1;
end