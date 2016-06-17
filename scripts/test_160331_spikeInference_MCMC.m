%% script - test spike inference MCMC with artifact removal
gcamp = 'f';
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
param.showFlag = 1;
[sig_mat_ra, a_mat_corrected] = removeArtifact( sig_mat, fName, param ); % removed artifact
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
plot([norm_nc(sig_mat(:,k),5) norm_nc(trend(:,k),5)]);axis tight
subplot(312);
plot([norm_nc(sig_mat_ra(:,k)-trend(:,k),5) norm_nc(bl(:,k),5)]);axis tight
subplot(313);
plot(df_mat(:,k),'r');axis tight
%% MCMC

% initialize params
disp([' ********* Spike Inference for GCaMP6',gcamp,' ********* ']);
k = 22; % 34, 10, 2, 4, 9 ,48
params.f = fs;          % time step size
params.prec = 5e-2; % prec specifies to what extent you want to discard the long slowly decaying tales of the ca response
params.B = 150;            % number of burn in samples (default 200)
params.Nsamples = 300;    % number of samples after burn in (default 500)
params.b = 0; % initializer for baseline
params.p = 2;                 % order of AR model (p == 1 or p == 2, default 1)
params.marg = 0;         % flag for marginalized sampler (default 0)
params.upd_gam = 1;      % flag for updating gamma (default 0) gamma = 1-dt/tau
if(strcmp( gcamp,'s'))
    params.A_lb = 0.115 * range(df_mat(:,k));
    params.tau_minmax = [0 1 3.5 15]; % GCAMP6S min and max raise time, min and max decay time
    [tau_rise,tau_decay] = deal(0.2 , 6); % GCAMP6S in sec
    params.TauStd = [0.05 0.5].*params.f; % standard deviation from time constants in samples!
elseif(strcmp( gcamp,'f')) 
    params.A_lb = 0.1 * range(df_mat(:,k));
    params.tau_minmax = [0 2 0.15 2]; % GCAMP6S min and max raise time, min and max decay time
    [tau_rise,tau_decay] = deal(0.1 , 1.5); % GCAMP6S in sec
    params.TauStd = [0.01 0.5].*params.f; % standard deviation from time constants in samples!
end
[g2,h] = tau_c2d(tau_rise,tau_decay,dt);
params.g = g2;
x = df_mat(:,k);
SAMPLES = cont_ca_sampler(x,params);  

C = make_mean_sample(SAMPLES,x);
S = mean(samples_cell2mat(SAMPLES.ss,size(C,2)));
P.b = mean(SAMPLES.Cb);
P.c1 = mean(SAMPLES.Cin);
P.neuron_sn = sqrt(mean(SAMPLES.sn2));
P.gn = mean(exp(-params.f./SAMPLES.g));
P.samples_mcmc = SAMPLES; % FN added, a useful parameter to have.

tau_ = mean(SAMPLES.g * dt)';
P.loglikeli = -norm(x' - C)^2;
disp('Time Constanst estimated :');
disp(['tau_rise is - ',num2str(tau_(1)),' [sec]']);
disp(['tau_decay is - ',num2str(tau_(2)),' [sec]']);
disp(['negtive log likelihood is : ',num2str(P.loglikeli)]);
%% Plot 

options.lineWidth = 1;
options.artifact = a_mat_corrected(:,k);
options.baseline = P.b;
F_raw = sig_mat(:,k);
F = x;
S(S < 0.1) = 0;
showInference(time, F_raw, F, C, S, options);
disp(['loglikeli of calcium trace for cell ',num2str(k),' is : ', num2str(P.loglikeli)]);
disp(['tau_rise is - ',num2str(tau_(1)),' [sec]']);
disp(['tau_decay is - ',num2str(tau_(2)),' [sec]']);
clear options;
%%
plot_continuous_samples(SAMPLES,x);
%% spike inference on complete dataset

P_batch = cell(size(df_mat,2),1);
params.print_flag = 1;
parfor j = 1:size(df_mat,2)
%     params.A_lb = 0.115 * range(df_mat(:,j));
    f = (df_mat(:,j));
    SAMPLES = cont_ca_sampler(f,params);  
    P_batch{j}.C = make_mean_sample(SAMPLES,f);
    P_batch{j}.S = mean(samples_cell2mat(SAMPLES.ss,size(C,2)));
    P_batch{j}.b = mean(SAMPLES.Cb);
    P_batch{j}.c1 = mean(SAMPLES.Cin);
    P_batch{j}.neuron_sn = sqrt(mean(SAMPLES.sn2));
    P_batch{j}.gn = mean(exp(-params.f./SAMPLES.g));
    P_batch{j}.samples_mcmc = SAMPLES; % FN added, a useful parameter to have.
    P_batch{j}.tau = mean(SAMPLES.g * dt)';
    P_batch{j}.loglikeli = -norm(f' - P_batch{j}.C)^2;
end
params.print_flag = 0;
%% Spike Thresholding

% simple senario  > 0.1
for trail = 1:size(df_mat,2)
    P_batch{trail}.S(P_batch{trail}.S < 0.1) = 0;
end
%% plot matrix

c_mat = cellfun(@(x) x.C,P_batch,'UniformOutput',false);
% matplot_examples = [1, 4, 9, 10, 11, 16, 18, 23, 26, 28, 29, 32, 38, 43, 45, 49];
matplot_examples = (1:size(sig_mat,2));
c_mat = cell2mat(c_mat)';
figure(9); 
subplot(121);imagesc(c_mat(:,matplot_examples)');
title('Calcium estimated matrix');ylabel('time [s]');
subplot(122);imagesc(df_mat(:,matplot_examples)');
title('F_{raw} matrix');ylabel('time [s]')

s_mat_norm = cellfun(@(x) z1(x.S),P_batch,'UniformOutput',false);
s_mat_norm = sparse(cell2mat(s_mat_norm)');
s_mat = cellfun(@(x) x.S,P_batch,'UniformOutput',false);
s_mat = sparse(cell2mat(s_mat)');

figure(10);
subplot(122);spyc(s_mat_norm(:,matplot_examples)','Parula',1);
daspect([170 2 1])
title('Raster plot');ylabel('cell');xlabel('time')
subplot(121);imagesc(c_mat(:,matplot_examples)');
colorbar;
daspect([170 2 1])
title('F_{raw} matrix');ylabel('time [s]')
%% plot single trace
ii = 1; % 34, 10, 2, 4, 9 ,48
options.lineWidth = 1;
F_raw = sig_mat(:,ii);
F = df_mat(:,ii);
C = P_batch{ii}.C;
S = P_batch{ii}.S;
S(S < 0.1) = 0;
options.baseline = P_batch{ii}.b;
if(ii ~= 1)
    options.artifact = a_mat_corrected(:,mod( ii,7 )+1);
else
    options.artifact = a_mat_corrected(:,1);
end

showInference(time, F_raw, F, C, S,options);
disp(['loglikeli of calcium trace for cell ',num2str(ii),' is : ', num2str(P_batch{ii}.loglikeli)]);
disp(['tau_rise is - ',num2str(P_batch{ii}.tau(1)),' [sec]']);
disp(['tau_decay is - ',num2str(P_batch{ii}.tau(2)),' [sec]']);
%%
plot_continuous_samples(P_batch{ii}.samples_mcmc,df_mat(:,ii));
%% plot comparison of response decay

temp = reshape(P_batch, [],4);
l = 4;
plotNum = 1;
for ii = (l:7:28)
    lw = 1;
    F_raw = sig_mat(:,ii);
    F = df_mat(:,ii);
    C = P_batch{ii}.C;    
    S = P_batch{ii}.S;
    S(S < 0.1) = 0;
    B = P_batch{ii}.b;
    if(ii ~= 1)
        arti = a_mat_corrected(:,mod( ii,7 )+1);
    else
        arti = a_mat_corrected(:,1);
    end

    subplot(4,1,plotNum);
    S(S == 0) = NaN;
    h1 = stem(time,S,'Color',[0.494117647409439 0.494117647409439 0.494117647409439]);
    h1.ShowBaseLine ='off';
    h1.BaseValue = B;
    hold on; 
    plot(time,F,'Color',[0.709215695858002 0.909803926944733 0.850980401039124]);
    title(['Spike inference with MCMC trial ',num2str(plotNum)],'FontSize',12);
    h2 = plot(time,C,'Color',[0.584313750267029 0.388235300779343 0.388235300779343],...
        'LineWidth',lw);
    h2.Color(4)=0.8;

    set(gca,'FontSize',13,'Box','off');
    xlabel('Time [sec]');
    ylabel('normalized C');
    if(any(arti))
        lw_fill = B * ones(size(arti,1),1);
        stim_logic = ones(size(arti,1),1);
        lw_fill(arti <  0.1)= B;
        stim_logic(arti <  0.5) = B;
        if(size(time,1) > size(time,2))
            h4 = fill( [time' fliplr(time')],  [stim_logic' lw_fill'], 'b');
        else 
            h4 = fill( [time fliplr(time)],  [stim_logic lw_fill], 'b');
        end
        alpha(.09);
        set(h4,'EdgeColor',[1 1 1]);
        set(h4,'EdgeAlpha',0);
        legend('F','est. C','est n_t','stimulus');
        set( legend,'Position',[0.855546467508258 0.432856164761928 0.10220768709019 0.197000685953818] );
        axis([time(1) time(end) (min([F(:) ;C(:)])-eps) max([(-1 * S(:));  F(:) ;C(:)])+eps]);
    else
        legend('F','est. C','est n_t');
        set( legend,'Position',[0.855546467508258 0.432856164761928 0.10220768709019 0.197000685953818] );
        axis('tight')
    end
    hold off
    plotNum = plotNum+1;
end
%% show decay graph

% for every cell (out of 7 cells)
trail_ind = reshape((1:28), [], 4);
for trail = 1:4
    ti_iter = trail_ind(:,trail);
    C = cell2mat(cellfun(@(x) x.C, P_batch(ti_iter),'UniformOutput',false));
    S = cell2mat(cellfun(@(x) x.S, P_batch(ti_iter),'UniformOutput',false));
    S(S < 0.1) = 0;
    S = logical(S);
    
    stim_logic = mean(a_mat_corrected,2);
    stim_logic(stim_logic <  0.1) = 0;
    stim_logic(stim_logic >=  0.1) = 1;
    a = diff(stim_logic);
    a(1) = 1; a(end) = 0;
    on_ind = find(abs(a));
    on_ind = reshape(on_ind,2,[]);
    C(~S) = NaN;
    for n = 1:size(on_ind,2)
        b = on_ind(:,n);
        C_on = C(:,(b(1):b(2)));
        C_box{trail}{n} = C_on(~isnan(C_on));
    end
end

trail_times = [{'06:51 AM'};{'07:00 AM'};{'07:10 AM'};{'07:27 AM'}];
for n = [0,2]
    figure;
    subplot(2,1,1)
    C_trail = {C_box{1+n,:}}';
    grp_size = cellfun(@numel,C_trail)';
    label = (1:numel( grp_size) ) ;
    grp = cell2mat( arrayfun(@(x,y) y * ones(1,x),grp_size, label,'UniformOutput',false) );
    boxplot(cell2mat(C_trail),grp)
    title(['Trail: ',num2str(n+1),' time: ',trail_times{n+1}])
    ylabel('$$\hat{C}$$','Interpreter','Latex');
    xlabel('Time ($$\Delta$$ t = 10 sec)','Interpreter','Latex');
    set(gca,'FontSize',11,'YLim',[0,0.3]);

    subplot(2,1,2)
    C_trail = {C_box{2+n,:}}';
    grp_size = cellfun(@numel,C_trail)';
    label = (1:numel( grp_size) ) ;
    grp = cell2mat( arrayfun(@(x,y) y * ones(1,x),grp_size, label,'UniformOutput',false) );
    boxplot(cell2mat(C_trail),grp)
    title(['Trail: ',num2str(n+2),' time: ',trail_times{n+2}])
    ylabel('$$\hat{C}$$','Interpreter','Latex');
    xlabel('Time ($$\Delta$$ t = 10 sec)','Interpreter','Latex');
    set(gca,'FontSize',11);
    set(gca,'FontSize',11,'YLim',[0,0.3]);
end

% trail_times_num_spikes = [{['06:51 AM',char(10),'num Spikes: 412']};{['07:00',char(10),'num Spikes: 575']};...
%     {['07:10 AM',char(10),'num Spikes: 97']};{['07:27 AM',char(10),'num Spikes: 0']}];
trail_times_num_spikes = [{'06:51 AM'};{'07:00 AM'};{'07:10 AM'};{'07:27 AM'}];
C_box_t = C_box';
C_ovr = cell2mat(C_box_t(:));
grp_size = sum(cellfun(@numel,C_box_t));
label = (1:numel( grp_size) ) ;
grp = cell2mat( arrayfun(@(x,y) y * ones(1,x),grp_size, label,'UniformOutput',false) );
figure
boxplot(C_ovr,grp)
title(['Overall estimated Ca^{2+} conc. at spikestimes vs. Trail'])
ylabel('$$\hat{C}$$','Interpreter','Latex');
xlabel('Trail','Interpreter','Latex');
set(gca,'FontSize',11,'YLim',[0,0.3],'XTickLabel',trail_times_num_spikes);
%%   Variables:
%   y:      raw fluorescence data (vector of length(T))
%   c:      denoised calcium concentration (Tx1 vector)
%   b:      baseline concentration (scalar)
%  c1:      initial concentration (scalar)
%   g:      discrete time constant(s) (scalar or 2x1 vector)
%  sn:      noise standard deviation (scalar)
%  sp:      spike vector (Tx1 vector)
%% params.g              discrete time constant(s) (estimated if not provided)
% params.sn             initializer for noise (estimated if not provided)
% params.b              initializer for baseline (estimated if not provided)
% params.c1             initializer for initial concentration (estimated if not provided)
% params.marg           flag for marginalized sampler (default 0)
% params.upd_gam        flag for updating gamma (default 0)
% params.gam_step       number of samples after which gamma is updated (default 50)
% params.std_move       standard deviation of shifting kernel (default 3*Dt)
% params.add_move       number of add moves per iteration (default T/100)
% params.init           initial sample 
% params.f              imaging rate (default 1)
% params.p              order of AR model (p == 1 or p == 2, default 1)
% params.defg           default discrete time constants in case constrained_foopsi cannot find stable estimates
% params.TauStd         standard deviation for time constants in continuous time (default [0.2,2])
