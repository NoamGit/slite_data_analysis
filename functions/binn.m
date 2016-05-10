function [ y ] = binn( X, binsize, method )
%binn( X, binsize, method ) takes vector X and biins its values according
%to the specified method and the binnsize
% method = @sum/ @mean/ @median etc..

if(size(X,1) == 1) 
    X = X';
end 
l = length(X);
numBins = floor( l/binsize );
y = reshape(X, binsize , numBins)';
y = method(y,2);
end

