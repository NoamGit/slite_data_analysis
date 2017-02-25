function [  ] = parseData_fn( fName, pathName )
%parseData_fn gets file path string and saves the processed data into a
%Parsing directory in the original location. This is very similar to the
%script parseData.m
% author: Noam Cohen, 2016

%% load data

pathTokens = textscan( pathName, '%s','delimiter',{' ','_','\','#','  ','Ret'});

if(~isempty(find(strcmp( 'gcamp6s',pathTokens{1} ),1)))
    gcamp = 's';
elseif(~isempty(find(strcmp( 'gcamp6f',pathTokens{1} ),1)))
    gcamp = 'f';
end

df = readtable([pathName fName]); % data frame
time = df.Time;
dt = double(time(2)-time(1));
fs = 1/dt;
off_period = 15*fs;
stim = df.Stim;
str = textscan(fName,'%s','delimiter','_');
full_field_s = [{'CRP'},{'WGN'}];
isnew_flag = isempty( findstr(pathName,'old') );
if(any(ismember(full_field_s,str{1})))
    sig_table = df(off_period:end,(3:end-1)); 
    dr = df(off_period:end,end); % dark region in the FOV
    full_field = 1;
else
    sig_table = df(off_period:end,(3:end)); 
end

num_sig = size(sig_table,2);
time = time(1:size(sig_table,1));
cellNames = [df.Properties.VariableNames];
cellNames = cellNames(3:end);
getCoord = @(txt) textscan(txt,'%f','delimiter','_');
getCoord_old = @(txt) textscan(txt,'%f','delimiter','-');
numInName = textscan(cellNames{1},'%s','delimiter','_');
cell_text = (cellfun(@(x) textscan(x,'%s','delimiter','_') ,cellNames)); % format 'Cell_cellIndx_locX_locY'
cellIndex = cell2mat(cellfun(@(c) str2double(c{2}),cell_text,'UniformOutput',false)); 
cellCenters = (cellfun(@(c)...
    [str2double(c{3}),str2double(c{4})],cell_text,'UniformOutput',false)); % [x ; y]

% artifact removal
param.showFlag = 0;
param.wavelets_levels = 3;
param.off_period = off_period;
param.showFlag = 0;
param.pathName = pathName;
[sig_table_ra, a_table_corrected] = removeArtifact( sig_table, fName, param ); % removed artifact
%% DF (data is after dark noise removal)
    
Fs_filt = 10;  % Sampling Frequency
N     = 5;     % Order
Fstop = 0.02;  % Stopband Frequency
Astop = 80;    % Stopband Attenuation (dB)
h = fdesign.lowpass('n,fst,ast', N, Fstop, Astop, Fs_filt);
Hd = design(h, 'cheby2', 'SystemObject', true);

parfor k = 1:num_sig
    
    x = sig_table_ra{:,k};
    
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
    Df_mat(:,k) = (x_detrend-bl(:,k))./abs(bl(:,k));
end

Df_table = array2table(single(Df_mat),'VariableNames',sig_table_ra.Properties.VariableNames);
%% MCMC

% initialize params
% disp([' ********* Spike Inference for GCaMP6',gcamp,' ********* ']);
k = 1; % 34, 10, 2, 4, 9 ,48
params.f = fs;          % time step size
params.prec = 5e-2; % prec specifies to what extent you want to discard the long slowly decaying tales of the ca response
params.B = 150;            % number of burn in samples (default 200)
params.Nsamples = 300;    % number of samples after burn in (default 500)
% params.B = 1;            % number of burn in samples (default 200)
% params.Nsamples = 1;   
params.b = 0; % initializer for baseline
params.p = 2;                 % order of AR model (p == 1 or p == 2, default 1)
params.marg = 0;         % flag for marginalized sampler (default 0)
params.upd_gam = 1;      % flag for updating gamma (default 0) gamma = 1-dt/tau
if(strcmp( gcamp,'s'))
    A_lb_init = 0.115;
    params.tau_minmax = [0 1 3.5 15]; % GCAMP6S min and max raise time, min and max decay time
    [tau_rise,tau_decay] = deal(0.2 , 6); % GCAMP6S in sec
    params.TauStd = [0.05 0.5].*params.f; % standard deviation from time constants in samples!
elseif(strcmp( gcamp,'f')) 
    A_lb_init = 0.1;
    params.tau_minmax = [0 2 0.15 2]; % GCAMP6S min and max raise time, min and max decay time
    [tau_rise,tau_decay] = deal(0.1 , 1.5); % GCAMP6S in sec
    params.TauStd = [0.01 0.5].*params.f; % standard deviation from time constants in samples!
end
[g2,h] = tau_c2d(tau_rise,tau_decay,dt);
params.g = g2;
% x = Df_mat(:,k);
P_batch = cell(size(Df_mat,2),1);
params.print_flag = 1;
parfor j = 1:size(Df_mat,2)
%     params.A_lb = 0.115 * range(df_mat(:,j));
    try
        params_iter = params;
        params_iter.A_lb = A_lb_init * range(Df_mat(:,k));
        f = (Df_mat(:,j));
        SAMPLES = cont_ca_sampler(f,params);  
        P_batch{j}.C = make_mean_sample(SAMPLES,f);
        P_batch{j}.S = mean(samples_cell2mat(SAMPLES.ss,size(P_batch{j}.C,2)));
        P_batch{j}.b = mean(SAMPLES.Cb);
        P_batch{j}.c1 = mean(SAMPLES.Cin);
        P_batch{j}.neuron_sn = sqrt(mean(SAMPLES.sn2));
        P_batch{j}.gn = mean(exp(-params.f./SAMPLES.g));
        P_batch{j}.samples_mcmc = SAMPLES; % FN added, a useful parameter to have.
        P_batch{j}.tau = mean(SAMPLES.g * dt)';
        P_batch{j}.loglikeli = -norm(f' - P_batch{j}.C)^2;
    catch exception
        disp(['inference problen in cell- ',num2str(j),' because:']); 
        disp(exception);
    end
end
params.print_flag = 0;
%% organize and save data to *.mat file

stim_list = [{'ORI'};{'BDN'};{'CNT'};{'CRP'};{'WGN'}];
stim_idx = (ismember(str{1},stim_list));
date_idx = [false; strcmp('data', pathTokens{1})];
date_idx = date_idx(1:end-1);
date_string = cell2mat(pathTokens{1}(date_idx));
Y = str2double(['20',date_string(1:2)]);M = str2double(date_string(3:4));D = str2double(date_string(5:6));

properties.stim_type = str{1}{stim_idx};
properties.date = datetime(Y,M,D);
properties.gcampType = gcamp;
properties.imaging_rate = fs;
properties.numCell = size(Df_mat,2);
properties.signalLenght = size(Df_mat,1);
properties.isnew = isnew_flag;
T = table;

for k = 1:properties.numCell
    if ~isempty(P_batch{k})
        frame.raw = sig_table{:,k};
        if(~isempty(a_table_corrected))
            frame.artifact = a_table_corrected{:,k};
        else 
            frame.artifact = 0;
        end
        frame.location = cellCenters{:,k};
        frame.Df = Df_mat(:,k);
        frame.C = P_batch{k}.C;
        frame.S = P_batch{k}.S;
        frame.tau = P_batch{k}.tau;
        frame.ll = P_batch{k}.loglikeli;
        frame.mcmc_samples = P_batch{k}.samples_mcmc;
        head_k = ['Cell_',num2str( cellIndex(k) )];
        T1 = table(frame, 'VariableName',{head_k} );
        T = [T T1];
    end
end

dataframe.properties = properties;
dataframe.cell_table = T;
txt = textscan( fName, '%s','delimiter',{'_sub','_',});
filename = cat(1, txt{1}(2:end-1), datestr(properties.date));
mkdir([pathName,'\Parsing\'])
cd([pathName,'\Parsing\'])
save( strjoin(filename,'_'), 'dataframe' );
end

