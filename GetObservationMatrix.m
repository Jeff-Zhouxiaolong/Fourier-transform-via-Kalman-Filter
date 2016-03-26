%function to generate observation matrix based on number of sample points
%much like the IDFT matrix but real, taking advantage of the fact that the
%input signal is real so IDFT has to produce real value
function [H] = GetObservationMatrix(Num_sample_points, A, H)
for j=1:Num_sample_points
    %put the matrix between 2->N
    counter = 1;
    for i=1:(Num_sample_points-1)
        theta = 2*pi*i*(j-1)/Num_sample_points;
        %put this in intermediate matrix between joining
        A(counter) = cos(theta);
        counter  = counter + 1;
        A(counter) = -sin(theta);
        counter = counter +1;
    end
    H(j,:) = [1 0 A];
    %correct location of N/2
%     H(j,Num_sample_points+1) = cos(pi*(j-1));
%     H(j,Num_sample_points+2) = 0;
    %actual matrix is [0 1 (2 till N/2-1) N/2 (2 till N/2 -1)]
end
H = H/Num_sample_points;
