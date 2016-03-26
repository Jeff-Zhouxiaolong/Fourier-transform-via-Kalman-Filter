%just checking the sub-sampling and sub-precision method
clc; clear all;
Num_sample_points_original = 16;
[A] = GetObservationMatrix(Num_sample_points_original);

Num_sample_points = 8;

index_location = zeros(1, 2*Num_sample_points);
index_location(1:2) = [1 2];

index = 2;
count = 1;
factor = 2; %division of new duration

while(index<2*Num_sample_points)
    index_location(index+1:index+2) = [count + 2*factor count+2*factor+1];
    count = count + 2*factor;
    index = index + 2;
end

display(A(:,index_location));

[B] = GetObservationMatrix(Num_sample_points);
display(B);

display(isequal(A(1:Num_sample_points_original/2, index_location), B));

A(1:Num_sample_points_original/2, index_location)== B