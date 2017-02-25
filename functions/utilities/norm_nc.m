function [ y ] = norm_nc( x, varargin )
% normalizes signal according to specified method
%   method 1 - according to abs(max) value
%   method 2 - according to max value
%   method 3 - according to min value
%   method 4 - according to std value
%   method 5 - rescale to [0,1] by X_new = ( x-Xmin )/( Xmax-Xmin )
%   method 6 - standartize (zero mean and 1 var) X_new = ( x-mean(x) )/( std(x) )


if nargin > 1
    
    if(numel(size(x)>1))
        x = x(:);
    end
    method = varargin{1};
    Xmax = max(x); 
    Xmin = min(x); 
    mu = mean(x); 
    sigma = std(x);
    
    switch method
        case 0 
            y = x;
        case 1 
            y = x./max(abs(x));
        case 2 
            y = x./Xmax;
        case 3 
            y = x./Xmin;
        case 4 
            y = x./sigma;
        case 5
            y = ( x-Xmin )./( Xmax-Xmin );
        case 6 
            y = ( x- mu )./( sigma );
        otherwise 
            error('Normalization method exception: method index is not compatible');
    end
else
    y = x./max(abs(x)); % default case
end

end