%% script - test spike inference constrained-foopsi with artifact removal

[fName, pathName] = uigetfile('D:\# Projects (Noam)\# SLITE\# DATA\*.xlsx','Select Cell''s Data Sheet');
%% load data

sheet_cell = xlsread([pathName fName]);
time = sheet_cell(:,1);
dt = double(time(2)-time(1));
fs = 1/dt;
off_period = 15*fs;
stim = sheet_cell(:,2);
str = textscan(fName,'%s','delimiter','_');
full_field_s = [{'CRP'},{'WGN'}];
if(any(ismember(full_field_s,str{1})))
    sig_mat = sheet_cell(off_period:end,(3:end-1)); 
    dr = sheet_cell(off_period:end,end); % dark region in the FOV
    full_field = 1;
else
    sig_mat = sheet_cell(off_period:end,(3:end)); 
end
num_sig = size(sig_mat,2);
time = time(1:size(sig_mat,1));

% artifact removal
param.showFlag = 0;
param.wavelets_levels = 3;
param.off_period = off_period;
sig_mat_ra = removeArtifact( sig_mat, fName, param ); % removed artifact
%% DF (data is after dark noise removal)
    
Fs_filt = 10;  % Sampling Frequency
N     = 5;     % Order
Fstop = 0.02;  % Stopband Frequency
Astop = 80;    % Stopband Attenuation (dB)
h = fdesign.lowpass('n,fst,ast', N, Fstop, Astop, Fs_filt);
Hd = design(h, 'cheby2', 'SystemObject', true);

for k = 1:num_sig
    
    x = sig_mat_ra(:,k);
    
    % step 1 - removing low frq trend with LPF
    trend(:,k) = filtfilt(Hd.SOSMatrix,Hd.ScaleValues,x);
    x_detrend = x - trend(:,k); % x1 is the detrended signal
    
    % step 2 -  cut 15 first sec from data and artifact
    x_detrend = x_detrend + mean(x); % avoiding deviding by 0 with bl
    
    % step 3 - Baseline estimation according lowest 15%
    Delta_t = 50 * fs; % sec * sr
    s_size = size(x_detrend,1);
    T = time(end);
    num_bin = floor( s_size / Delta_t );
    s_mat_binned = reshape(x_detrend(1:num_bin*Delta_t,:),Delta_t, num_bin, []);
    t_mat_binned = reshape(time(1:num_bin*Delta_t,:),Delta_t, num_bin, []);
    perc = .20; % the percent that the hist values under are probabely bsl
    quantile_vec = quantile(s_mat_binned, perc, 1);
    indx_mat = bsxfun(@(a,b) a < b,s_mat_binned,quantile_vec);
    Y_val = s_mat_binned(indx_mat); X_val = t_mat_binned(indx_mat);
    p = polyfit( X_val, Y_val,1);
    bl(:,k) = polyval(p,time);
    df_mat(:,k) = (x_detrend-bl(:,k))./abs(bl(:,k));
end
%% plotting test

k = 4;
figure(6)
subplot(311);
plot([norm_nc(sig_mat_ra(:,k),5) norm_nc(trend(:,k),5)]);axis tight
subplot(312);
plot([norm_nc(sig_mat_ra(:,k)-trend(:,k),5) norm_nc(bl(:,k),5)]);axis tight
subplot(313);
plot(df_mat(:,k),'r');axis tight
%% c-foopsi

k = 4; % 34, 10, 2, 4, 9 ,48 ||     noise -35
x = df_mat(:,k); % raw fluorescence data (vector of length(T))

% initialize params
%   b:      baseline concentration (scalar)
%  c1:      initial concentration (scalar)
%   g:      
%  sn:      noise standard deviation (scalar)
%  sp:      spike vector (Tx1 vector)
options.dt = 0.1;
% g = exp(-options.dt./9); % discrete time constant(s) (scalar or 2x1 vector)
[tau_rise,tau_decay] = deal(0.2 , 9); % in sec
[g,h] = tau_c2d(tau_rise,tau_decay,dt);
options.p = 2; % order of AR model
sn = 1; % noise standard deviation (scalar)
options.method = 'spgl1';
options.resparse = 4;

[cc,cb,c1,gn,sn,spk] = constrained_foopsi(x,[],[],[],[],options);
                            
gd = max(roots([1,-gn(:)']));
gd_vec = gd.^((0:size(y,1)-1));
C = mean(bsxfun(@plus,cc',cb + c1*gd_vec),1);
S = mean(spk');
% C = full(cc(:)' + cb + c1*gd_vec);
C = cc(:,end)' + cb(end) + c1(end)*gd_vec;
S = spk(:,end)';
P.b = cb;
P.c1 = c1;           
P.neuron_sn = sn;
P.gn = gn;
tau_ = tau_d2c(gn,dt);
disp('Time Constanst estimated :');
disp(['tau_rise is - ',num2str(tau_(1)),' [sec]']);
disp(['tau_decay is - ',num2str(tau_(2)),' [sec]']);
%% Plot 

Pl.n    = S; 
Pl.n(Pl.n==0) = NaN; % true spike train (0's are NaN's so they don't plot)
Pl.lw   = 1;                    % line width

figure(8);
subplot(311);
h5 = plot(time,sig_mat(:,k),'LineWidth',Pl.lw);
axis tight;
title('F_{raw} (with arti)','FontSize',12);
set(gca,'FontSize',13,'Box','off');

subplot(312);
h3 = plot(time,(x),'k','LineWidth',Pl.lw);
title('dF_{processed}','FontSize',12);
h3.Color(4)=1;
hold on;
h4 = plot(time,C,'Color',[1 0.600000023841858 0.7843137383461]...
    ,'LineStyle','--','LineWidth',Pl.lw+0.5);
h4.Color(4)=0.8;
set(gca,'FontSize',14,'Box','off');
ylabel('\Delta{F}/F');
axis('tight')
hold off;

subplot(313);
plot(time,norm_nc(x,0),'Color',[0.709215695858002 0.909803926944733 0.850980401039124]);
title('Spike inference with constrained foopsi','FontSize',12);
hold on; 
h2 = plot(time,norm_nc(C,0),'Color',[0.584313750267029 0.388235300779343 0.388235300779343],...
    'LineWidth',Pl.lw);
h2.Color(4)=0.5;

h1 = stem(time,-1*Pl.n,'Color',[0.494117647409439 0.494117647409439 0.494117647409439]);
h1.ShowBaseLine ='off';
h1.Color(4)=0.01;
set(gca,'FontSize',13,'Box','off');
xlabel('Time [sec]');
ylabel('normalized C');
legend('F','est. C','est n_t');
set( legend,'Position',[0.855546467508258 0.432856164761928 0.10220768709019 0.197000685953818] );
axis('tight')
hold off

%   Variables:
%   y:      raw fluorescence data (vector of length(T))
%   c:      denoised calcium concentration (Tx1 vector)
%   b:      baseline concentration (scalar)
%  c1:      initial concentration (scalar)
%   g:      discrete time constant(s) (scalar or 2x1 vector)
%  sn:      noise standard deviation (scalar)
%  sp:      spike vector (Tx1 vector)

%   USAGE:
%   [c,b,c1,g,sn,sp] = constrained_foopsi(y,b,c1,g,sn,OPTIONS)
%   The parameters b,cin,g,sn can be given or else are estimated from the data

%   OPTIONS: (stuct for specifying options)
%         p: order for AR model, used when g is not given (default 2)
%    method: methods for performing spike inference
%   available methods: 'dual' uses dual ascent
%                       'cvx' uses the cvx package available from cvxr.com (default)
%                      'lars' uses the least regression algorithm 
%                     'spgl1' uses the spgl1 package available from
%                     math.ucdavis.edu/~mpf/spgl1/  (usually fastest)
%   bas_nonneg:   flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y)
%   noise_range:  frequency range over which the noise power is estimated. Default [Fs/4,Fs/2]
%   noise_method: method to average the PSD in order to obtain a robust noise level estimate
%   lags:         number of extra autocovariance lags to be considered when estimating the time constants
%   resparse:     number of times that the solution is resparsened (default 0). Currently available only with methods 'cvx', 'spgl'
%   fudge_factor: scaling constant to reduce bias in the time constant estimation (default 1 - no scaling)
