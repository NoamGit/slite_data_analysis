function [ delF, time ] = DeltaF( sig, span ,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% check data

% cases
if nargin < 3
    time = sig(:,1);
    delF = (sig(:,2))./(max(sig(:,2)) - min(sig(:,2)));
    delF = smooth(delF, span);

else
    backgr = varargin{1};
    time = sig(:,1);
    delF = (sig(:,2) - backgr(:,2))./(max(sig(:,2)) - min(sig(:,2)));
    delF = smooth(delF, span);

% different methods
if nargin > 3
    meth_IND = find(strcmp({varargin{2:end}},'bg_method'));
        switch varargin{meth_IND+2}
            case 'poly2'
                bg_delta_COEFF = polyfit( time ,delF, 2);
                bg_delta = polyval(bg_delta_COEFF, time);
                delF = delF - bg_delta + abs( min(delF - bg_delta) ) ; 
            otherwise
        end
end
end

