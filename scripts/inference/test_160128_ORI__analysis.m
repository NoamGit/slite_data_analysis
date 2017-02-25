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

k = 5;
subplot(311);
plot([norm_nc(sig_mat(15*fs:end,k),5) norm_nc(trend(15*fs:end,k),5)]);axis tight
subplot(312);
plot([norm_nc(sig_mat(15*fs:end,k)-trend(15*fs:end,k),5) norm_nc(bl(:,k),5)]);axis tight
subplot(313);
plot(df_mat(:,k),'r');axis tight
%% Artifact removal
% function input order(signal, arti, on_time, filt_bs, safety_offset, filt_thres ,fs, maxl_phase)

[ sm_proc, a_proc ] = ArtiRemOrientation(df_mat, arti_mat(15*fs-1:end,:), 5, 11, 30, 1e3, 10, 16);
figure(7); imagesc(sm_proc(:,[1:8,10:11,13:end])')
%% plotting after artifact removal 

k = 35;
figure(2);
subplot(411);
h1 = plot(t_trim , df_mat(:,k),'k');h1.Color(4)=1;axis tight;
hold on;
h2 = plot(t_trim , a_proc(:,k),'y');h2.Color(4)=0.4;axis tight;
hold off
legend('x','artifact')

subplot(412);
h3 = plot(t_trim, sm_proc(:,k));axis tight;h3.Color(4)=0.7;axis tight;
hold on;
h2 = plot(t_trim , a_proc(:,k),'y');h2.Color(4)=0.4;axis tight;
hold off
legend('x - artifact','artifact')

subplot(413);
h1 = plot(t_trim , df_mat(:,10),'k');h1.Color(4)=1;axis tight;
hold on;
h2 = plot(t_trim , a_proc(:,10),'y');h2.Color(4)=0.4;axis tight;
hold off

subplot(414);
h3 = plot(t_trim, sm_proc(:,10));axis tight;h3.Color(4)=0.7;axis tight;
hold on;
h2 = plot(t_trim , a_proc(:,10),'y');h2.Color(4)=0.4;axis tight;
hold off
xlabel('Time [sec]')
ylabel('DF')
%% for conference
% hand picking good indexes to show
matplot_examples = [1, 4, 9, 10, 11, 16, 18, 23, 26, 28, 29, 32, 38, 43, 45, 49];
figure(8); imagesc(sm_proc(:,matplot_examples)');
ylabel('time [s]')

% hand picking examples for full DF traces
df_examples = [2, 3, 4 ,10, 14];
figure(10);
subplot(511);plot(t_trim./60, sm_proc(:,matplot_examples(2))');axis('tight');
subplot(512);plot(t_trim./60, sm_proc(:,matplot_examples(14))');axis('tight');
subplot(513);plot(t_trim./60, sm_proc(:,matplot_examples(4))');axis('tight');
subplot(514);plot(t_trim./60, sm_proc(:,matplot_examples(10))');axis('tight');
subplot(515);plot(t_trim./60, sm_proc(:,matplot_examples(3))');axis('tight');
xlabel('Time [min]')
%% Kalman filtering after artifact removal (2DO)
% yv - obs; y - signal;  x - the dynamic state prcess

% Next, implement the filter recursions in a FOR loop:
P=B*Q*B';         % Initial error covariance
x=zeros(3,1);     % Initial condition on the state
ye = zeros(length(t),1);
ycov = zeros(length(t),1); 
errcov = zeros(length(t),1); 

for i=1:length(t)
  % Measurement update
  Mn = P*C'/(C*P*C'+R);
  x = x + Mn*(yv(i)-C*x);  % x[n|n]
  P = (eye(3)-Mn*C)*P;     % P[n|n]

  ye(i) = C*x;
  errcov(i) = C*P*C';

  % Time update
  x = A*x + B*u(i);        % x[n+1|n]
  P = A*P*A' + B*Q*B';     % P[n+1|n] 
end
%% engaging artifact in plots

artifact = mean( a_proc,2);
arti_std = stdfilt(artifact,ones(5,1));
plot(t_trim, norm_nc(artifact), t_trim,norm_nc(arti_std)>0.04);axis tight
arti_logical = norm_nc(arti_std) > 0.04;
