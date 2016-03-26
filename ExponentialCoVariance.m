%function return dimensionxdimension exponential co-variance matrix based on input parameters
%ini_std   -> standard deviation of all random variables, could be a vector
%dimension -> dimension of square matrix
%psi       -> variance of non-spatial error, measurement error or stochastic 
%             source error associated with each loation
%scale     -> determines the scale of exponential based on euclidean distance
%phi       -> range of spatial dependence

function [Covariance] = ExponentialCoVariance(ini_std,dimension, psi, scale, phi)
%put the main co-variances in the diagonal location
v = ones(1,dimension).*ini_std.*ini_std;
Covariance = diag(v) + psi*eye(dimension);
%calculate the rest of the locations, following exponential co-variance formula
scale = ini_std*ini_std;

for i=1:dimension
    for j=1:dimension
        if i==j,
            continue;
        end
        Covariance(i,j) = scale*exp(-abs(i-j)/phi);
    end    
end

