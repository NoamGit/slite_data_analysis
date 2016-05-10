function [ y, x_trend ] = detrend( x, param )
% detrends signal according to parameter

    if(nargin < 2)
        Fs_filt = 10;  % Sampling Frequency
        N     = 5;     % Order
        Fstop = 0.035;  % Stopband Frequency
        Astop = 80;    % Stopband Attenuation (dB)
    else
        if isfield(param,'N'); N = param.N;end
        if isfield(param,'Fstop'); Fstop = param.Fstop; end
        if isfield(param,'Astop'); Astop = param.Astop; end
        if isfield(param,'Fs_filt'); Fs_filt = param.Fs_filt; end
    end

    h = fdesign.lowpass('n,fst,ast', N, Fstop, Astop, Fs_filt);
    Hd = design(h, 'cheby2', 'SystemObject', true);
    x_trend = filtfilt(Hd.SOSMatrix,Hd.ScaleValues,x);
    y = x - x_trend; % x1 is the detrended signal

end

