function [ S_tau, fs, fft_tau ] = ZOHSpectralEst( in, dt )
% perform spectral analysis of point proccesses using ZOH interpulation
% dummy data :
% t = (0:0.01:5);
% sp = [ 0.34 0.6 2.44 3.21 3.22 4 4.45 4.6 4.88];
% tau = double(ismember(t,sp)); v = [1 1 3 1 2 1 2 1 ];
% tau(logical(tau)) = v;
%% constants 

NW = 1.5;
precent = 0.5;
%% creates the Tau-process from sp (can have values > 1)
% binarize if needed
sp = in;
if(any(arrayfun(@(x) rem(x,1), nonzeros(sp))))
    sp(sp < quantile(sp,precent)) = 0; sp = ceil(sp);
end

% build new signal with new resolution accordint to max sp
T = (length(sp)-1) * dt;
m = 1;
sp_r = zeros(m * length(sp),1); 
sp_r(m * find(sp)) = nonzeros(sp);
if numel(unique(sp))>2
    for k = find(sp_r > 0)'
        extrap_factor = sp_r(k);
        sp_r(k-extrap_factor+1:k) = 1;
    end
end
t_r = (0:dt/m:(T+dt)-dt/m);
sp_times = t_r(find(sp_r));
tau_val = [sp_times(1) diff(sp_times) t_r(end)-sp_times(end)];
tau = sp_r; indx = [1 ;find(sp_r)]; tau(indx(1:end)) = tau_val;
% subplot(211);stem(t,sp);subplot(212);stem(t_r,tau);axis tight; % sanity check
%% ZOH interpulation is equivilant to interp1 

t_nz = t_r(logical(tau));
tau_nz = tau(logical(tau));
tau_zoh = interp1(t_nz', tau_nz, t_r,'previous');
tau_zoh(isnan(tau_zoh)) = tau_nz(end); % some hack...FIXME
% subplot(211);bar(t_r,tau_zoh);axis tight; % sanity check
%% estiamte frequency response

fs = 1/(dt/m);
Ns = length(sp);
df = fs/Ns;
f = ( -fs/2:df:(fs/2-df) )';
% indx = f >=0; f = f(indx);
fft_tau =fft(tau_zoh);
% S_tau = S_tau(indx);

%  with MTM
S_tau = pmtm(tau_zoh, NW,Ns,fs);
% subplot(211);Hpsd = dspdata.psd(S_tau,'Fs',fs); plot(Hpsd)
% subplot(212);bar(t_r,tau_zoh);axis tight; % sanity check
end

