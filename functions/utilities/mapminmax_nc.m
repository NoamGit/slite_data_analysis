function [ y, fun ] = mapminmax_nc( x, ymin,ymax )
%mapminmax_nc( x, ymin,ymax ) Summary of this function goes here
%   Detailed explanation goes here
    xmax = max(x(:));
    xmin = min(x(:));
    if(xmax == xmin)
        y = ymin .*ones(length(x),1);
    else
        y = (ymax-ymin)*(x-xmin)/(xmax-xmin) + ymin;
        fun = @(X) (ymax-ymin)*(X-xmin)/(xmax-xmin) + ymin;
    end
end

