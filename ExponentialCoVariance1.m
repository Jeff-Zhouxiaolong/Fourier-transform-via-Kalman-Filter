%function return dimensionxdimension exponential co-variance matrix based on input parameters
%ini_std   -> standard deviation of all random variables, could be a vector
%dimension -> dimension of square matrix
%psi       -> variance of non-spatial error, measurement error or stochastic 
%             source error associated with each loation
%scale     -> determines the scale of exponential based on euclidean distance
%phi       -> range of spatial dependence

function [Covariance] = ExponentialCoVariance1(ini_std,dimension, psi, scale, phi)
%put the main co-variances in the diagonal location
v = ones(1,dimension).*ini_std.*ini_std;
Covariance = diag(v) + psi*eye(dimension);
%calculate the rest of the locations, following exponential co-variance formula
for i=0:dimension-1
    for j=0:dimension-1
        if i==j,
            continue;
        end
        k = i;
        l = j;
        %co-variance of conjugate thingy
%         if (i < dimension/2 +1) && (j>dimension/2 +1),
%             distance = abs(i - (j-(dimension/2)));
%         elseif (i > dimension/2 + 1) && (j<dimension/2 +1),
%             distance = abs((i-dimension/2) - j);
%         else
%             distance = abs(i-j);
%         end
        if(k > dimension/2 +1),
            k = k - dimension/2;
        end
        if(l > dimension/2 +1),
            l = l - dimension/2;
        end
        distance = abs(k - l);
        Covariance(i+1,j+1) = scale*exp(-distance/phi);
    end    
end

