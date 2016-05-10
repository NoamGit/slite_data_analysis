%% script - test spike inference from real data

[fName, pathName] = uigetfile('D:\# Projects (Noam)\# SLITE\# DATA\*.xlsx','Select Cell''s Data Sheet');
%% load data

sheet = xlsread([pathName fName]);
time = sheet(:,1);
dt = double(time(2)-time(1));
fs = 1/dt;
stim = sheet(:,2);
sig_mat = sheet(:,(3:end)); 
num_sig = size(sig_mat,2);
%% load artifact

[fName, pathName] = uigetfile('D:\# Projects (Noam)\# SLITE\# DATA\*.xlsx','Select Cell''s Data Sheet');
sheet = xlsread([pathName fName]);
arti_mat = sheet(:,(3:end)); 
%% DF (data is after dark noise removal)
    
Fs_filt = 10;  % Sampling Frequency
N     = 5;     % Order
Fstop = 0.02;  % Stopband Frequency
Astop = 80;    % Stopband Attenuation (dB)
h = fdesign.lowpass('n,fst,ast', N, Fstop, Astop, Fs_filt);
Hd = design(h, 'cheby2', 'SystemObject', true);

for k = 1:num_sig
    
    x = sig_mat(:,k);
    
    % step 1 - removing low frq trend with LPF
    trend(:,k) = filtfilt(Hd.SOSMatrix,Hd.ScaleValues,x);
    x_detrend = x - trend(:,k); % x1 is the detrended signal
    
    % step 2 -  cut 15 first sec from data and artifact
    x_detrend = x_detrend(15*fs:end);
    x_detrend = x_detrend + mean(x); % avoiding deviding by 0 with bl
    t_trim = time(15*fs:end);
    
    % step 3 - Baseline estimation according lowest 15%
    Delta_t = 50 * fs; % sec * sr
    s_size = size(x_detrend,1);
    T = t_trim(end);
    num_bin = floor( s_size / Delta_t );
    s_mat_binned = reshape(x_detrend(1:num_bin*Delta_t,:),Delta_t, num_bin, []);
    t_mat_binned = reshape(t_trim(1:num_bin*Delta_t,:),Delta_t, num_bin, []);
    perc = .15; % the percent that the hist values under are probabely bsl
    quantile_vec = quantile(s_mat_binned, perc, 1);
    indx_mat = bsxfun(@(a,b) a < b,s_mat_binned,quantile_vec);
    Y_val = s_mat_binned(indx_mat); X_val = t_mat_binned(indx_mat);
    p = polyfit( X_val, Y_val,1);
    bl(:,k) = polyval(p,t_trim);
    df_mat(:,k) = (x_detrend-bl(:,k))./abs(bl(:,k));
end
%% plotting test

k = 4;
figure(6)
subplot(311);
plot([norm_nc(sig_mat(15*fs:end,k),5) norm_nc(trend(15*fs:end,k),5)]);axis tight
subplot(312);
plot([norm_nc(sig_mat(15*fs:end,k)-trend(15*fs:end,k),5) norm_nc(bl(:,k),5)]);axis tight
subplot(313);
plot(df_mat(:,k),'r');axis tight
%% Artifact removal
% function input order(signal, arti, on_time, filt_bs, safety_offset, filt_thres ,fs, maxl_phase)

[ sm_proc, a_proc ] = ArtiRemOrientation(df_mat, arti_mat(15*fs:end,:), 5, 11, 30, 1e3, 10, 16);
sm_proc_norm = arrayfun(@(x) norm_nc(sm_proc(:,x)), (1:num_sig),'UniformOutput',false);
figure(7); imagesc(cell2mat(sm_proc_norm)')
%%

k = 10;
figure(6)
subplot(311);
plot([norm_nc(sig_mat(15*fs:end,k),5) norm_nc(trend(15*fs:end,k),5)]);axis tight
subplot(312);
plot([norm_nc(sig_mat(15*fs:end,k)-trend(15*fs:end,k),5) norm_nc(bl(:,k),5)]);axis tight
subplot(313);
plot(sm_proc(:,k),'r');axis tight
%% fast-oopsi

k = 34; % 34, 10, 2, 4, 9 ,48
x = sm_proc(:,k);

% set simulation metadata
V.F = norm_nc(x,5);
V.dt            = 0.1;       % time step size
V.fast_iter_max = 5;        % whether to plot with each iteration
V.fast_plot     = 1;
V.save          = 0;        % whether to save results

% initialize params
P.a     = max(x);                                           % scale
P.b     = 0.1;                                              % bias (baseline median)
tau     = 12;                                               % decay time constant for each cell
P.gam   = 1-V.dt./tau;                                      % set gam
P.lam   = 0.1;                                              % rate [mean rate]
P.sig   = median(abs(x - median(x)))/1.482;                 % F noise std
% P.sig   = 0.01;                % F noise std ADAPT

[n_best, P_best,~ , C] = fast_oopsi(V.F,V,P);
tau_best = -V.dt/(P_best.gam-1);

% gain factor determination (gfd) according to exp fit
kappa =60;
eta = 9;
G = (1/max(n_best)) *...
    (1-(1./(1+exp((-sum((z1(C)-z1(x)).^2)+kappa)./eta ))));
n_best_scaled = n_best * G;

% 3-level treshold
T=length(V.F);
if(P.sig > 0.03)
    tresh = quantile(n_best_scaled, 0.9925);
elseif(P.sig >= 0.025)
    tresh = quantile(n_best_scaled, 0.9925);
else
%     tresh = quantile(I{1}.n/max(I{1}.n), 0.999);
    tresh = 0.70;
end

%% Plot 

Pl.n    = double(n_best_scaled >= tresh ); 
Pl.n(Pl.n==0)=NaN; % true spike train (0's are NaN's so they don't plot)
Pl.lw   = 1;                    % line width
figure(8);
subplot(211);
plot(t_trim,(x),'Color',[0.494117647409439 0.494117647409439 0.494117647409439],'LineWidth',Pl.lw)
set(gca,'FontSize',14,'Box','off');
xlabel('Time [sec]');
ylabel('\Delta{F}/F');
axis('tight')
subplot(212);
plot(t_trim,z1(x),'Color',[0.709215695858002 0.909803926944733 0.850980401039124]);
hold on; plot(t_trim,z1(C),'Color',[0.584313750267029 0.388235300779343 0.388235300779343],...
    'LineWidth',Pl.lw)
stem(t_trim,Pl.n,'Color',[0.494117647409439 0.494117647409439 0.494117647409439]);
hold off
set(gca,'FontSize',14,'Box','off');
xlabel('Time [sec]');
ylabel('normelized C');
legend('F','est. C','est n_t');
set(legend,'Location','northwestoutside');
axis('tight')

