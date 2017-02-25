function [ y ] = binn( X, binsize, method )
%binn( X, binsize, method ) takes vector X and bins its values according
%to the specified method and the binnsize
% method = @sum/ @mean/ @median etc..

if(size(X,1) == 1) 
    X = reshape(X,[],1);
end 
l = length(X);
reminder = rem(l,binsize);
if(reminder)
    X = padarray(X,[reminder 0],'replicate','post');
end

numBins = length(X)/binsize;
X = reshape(X, binsize , numBins).';
y = method(X,2);
if size(y,2) > 1
    y= method(X,[],2);
end

end

