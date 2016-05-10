% Test artifact removal with spectral filtering
% requirements - 
%           in fullfield stimulation we need the last ROI to be an "empty" ROI
%           in spatiotemporal stimulation we need an artifact

[fName, pathName] = uigetfile('C:\Users\noambox\Dropbox\# MASTER\# SLITE\Data Analysis\scripts\artifact removal\*.xlsx','Select Cell''s Data Sheet');
%% load data

sheet_cell = xlsread([pathName fName]);
time = sheet_cell(:,1);
dt = double(time(2)-time(1));
fs = 1/dt;
off_period = 15*fs+1;
stim = sheet_cell(:,2);
str = textscan(fName,'%s','delimiter','_');
full_field_s = [{'CRP'},{'WGN'}];
full_field = 0;
if(any(ismember(full_field_s,str{1})))
    sig_mat = sheet_cell(off_period:end,(3:end-1)); 
    dr = sheet_cell(off_period:end,end); % dark region in the FOV
    full_field = 1;
else
    sig_mat = sheet_cell(off_period:end,(3:end)); 
end
num_sig = size(sig_mat,2);
time = time(1:size(sig_mat,1));
%% artifact removal

%*********************************************%
%**************** CONSTANTS ******************%
%*********************************************%
lev = 3;
        
% plot raw data
s = size(sig_mat,2);
figure(1);
for z = 1:s
    subplot(s,1,z)
    plot(time, sig_mat(:,z));
    ylabel('raw F data'); xlabel('time [sec]');
    axis tight
end

if(full_field)
    %% temporal artifact removal
    % with full field stimulus we only take the artifact from a dark region, traslate it to have zero
    % baseline and subtract it from all signals. 
    
    %*********************************************%
    %**************** CONSTANTS ******************%
    %*********************************************%
    param_f.Fs_filt = 10;  % Sampling Frequency
    param_f.N     = 5;     % Order
    param_f.Fstop = 0.035;  % Stopband Frequency
    param_f.Astop = 80;    % Stopband Attenuation (dB)
    win_size = 200;
    
    % detrend dark region
%     h = fdesign.lowpass('n,fst,ast', N, Fstop, Astop, Fs_filt);
%     Hd = design(h, 'cheby2', 'SystemObject', true);
%     dr_trend = filtfilt(Hd.SOSMatrix,Hd.ScaleValues,dr);
%     dr_detrend = dr - dr_trend; % x1 is the detrended signal
    dr_detrend = detrend( dr, param_f );
    dr_detrend = dr_detrend + abs(min(dr_detrend));
    
    % plot
    figure(2)
    subplot(311);plot(dr); axis tight
    subplot(312);plot(dr_trend); axis tight
    subplot(313);plot(dr_detrend); axis tight

    % Remove artifact and plot results
    param.win_size = win_size;
    [ sig_mat_t, a_mat ] = ArtiRemFullField(sig_mat , dr_detrend, param);
%     sig_mat_t = bsxfun(@minus, sig_mat, dr_detrend);
    
    % plot
    figure(2);
    range = (1:size(sig_mat_t,1)); % (2e3:3.5e3)
    k = 3;
    for z = 1:k
        subplot(k,1,z)
        h1 = plot(sig_mat(range,z));
        hold on; plot(sig_mat_t(range,z)); h1.Color(4) = 0.6;
        h2 = plot(-1 * a_mat(range,z) + min(sig_mat_t(range,z)),'k'); h2.Color(4) = 0.6;
        ylabel('F_t'); xlabel('time [sec]');
        legend('s','s - a',' - a');
        hold off;
        axis tight
    end
else
    %% spatio-temporal artifact removal
    if(any(ismember(full_field_s,'BDN'))) % case where we have binary dense noise
        sig_mat_t = wden(sig_mat,'sqtwolog','s','sln',lev,'sym8');
        sig_mat_t = reshape(sig_mat_t,size(sig_mat,1),[]);
        
        figure(3)
        range = (1:size(sig_mat_t,1)); % (2e3:3.5e3)
        for z = 1:4
            subplot(4,1,z)
            h1 = plot(sig_mat(range,z));
            hold on; plot(sig_mat_t(range,z)); h1.Color(4) = 0.5;
%             h2 = plot(-1 * a_mat(range,z) + min(sig_mat_t(range,z)),'k'); h2.Color(4) = 0.6;
            ylabel('F_t'); xlabel('time [sec]');
            legend('s','s filt');
            hold off;
            axis tight
        end
    else if(any(ismember(str{1},['CNT','ORI'])))
            
            %*********************************************%
            %**************** CONSTANTS ******************%
            %*********************************************%
            param_f.Fs_filt = 10;  % Sampling Frequency
            param_f.N     = 5;     % Order
            param_f.Fstop = 0.001;  % Stopband Frequency
            param_f.Astop = 80;    % Stopband Attenuation (dB)
            
            % load artifact
            disp('stimulus is spatially changing please load artifact');
            [fName, pathName] = uigetfile('D:\# Projects (Noam)\# SLITE\# DATA\*.xlsx','Select Cell''s Data Sheet');
            sheet_art = xlsread([pathName fName]);
            arti_mat = sheet_art(off_period:end,(3:end)); 
            
            % detrend art
%             h = fdesign.lowpass('n,fst,ast', N, Fstop, Astop, Fs_filt);
%             Hd = design(h, 'cheby2', 'SystemObject', true);
%             dr_trend = filtfilt(Hd.SOSMatrix,Hd.ScaleValues,arti_mat);
            arti_mat_d = arti_mat - dr_trend; % x1 is the detrended signal
            arti_mat_d = detrend(arti_mat , param_f);
 
            % plot
            figure(4)
            k = 4;
            subplot(311);plot(arti_mat(:,k)); axis tight
            subplot(312);plot(dr_trend(:,k)); axis tight
            subplot(313);plot(arti_mat(:,k) - dr_trend(:,k)); axis tight
            
            % plot 
            figure(5)
            for z = 1:4
                subplot(4,1,z)
                plot([ z1(sig_mat(:,z)) ]); hold on;
                h1 = plot(z1(arti_mat_d(:,z))); h1.Color(4) = 0.5;
                hold off
                ylabel('norm signal'); xlabel('time [sec]');
                legend('signal', 'artifact');
                axis tight
            end
            
            if(any(ismember(str{1},['CNT'])))
                %**************** CONSTANTS ******************%
                param.safty = 5;
                param.win_size = 25 * fs;
                [ sig_mat_t, a_mat ] = ArtiRemContrast(sig_mat, arti_mat_d, param);
            else
                %**************** CONSTANTS ******************%
                param.on_time = 5;
                param.std_filt_binsize = 11;
                param.safty = 30;
                param.filt_thres = 1e3;
                param.fs = 10;
                param.maxl_phase = 16; 
                [ sig_mat_t, a_mat ] = ArtiRemOrientation(sig_mat, arti_mat_d, param);
            end
           
            figure(6)
            range = (1:size(sig_mat_t,1)); % (2e3:3.5e3)
            for z = 1:4
                subplot(4,1,z)
                h1 = plot(sig_mat(range,z));
                hold on; plot(sig_mat_t(range,z)); h1.Color(4) = 0.5;
                h2 = plot(-1 * a_mat(range,z) + min(sig_mat_t(range,z)),'k'); h2.Color(4) = 0.6;
                ylabel('F_t'); xlabel('time [sec]');
                legend('s','s - a',' - a');
                hold off;
                axis tight
            end
        end
    end
end


% %%
% % extract noise variance
% fs = 10; % hz
% chirp_duration = 30; %sec
% y = sig_mat(:,z);
% noise = y(1:14*fs);
% s = y(14*fs:end);
% % sigma = real(fft(noise)).^2;
% noise_sigma = var(noise);
% reps = floor(length(y(15*fs:end))/(chirp_duration*fs));
% arti_sigma = var(y(15*fs:15*fs + chirp_duration*fs*reps));
% PSF = reshape(y(15*fs:15*fs - 1 + chirp_duration*fs*reps),chirp_duration*fs,[]);
% PSF = bsxfun(@minus, PSF ,mean(PSF,1));
% for k = 2:size(PSF,2)
%     PSF(:,k) = correctPhase(PSF(:,1),PSF(:,k),100);
% end
% PSF_i = PSF(:,2);
% PSF_i = PSF(:,3);
% 
% % wiener filter
% wnr = deconvwnr(s, PSF_i,noise_sigma/arti_sigma);
% plot(conv(wnr,flipud(PSF(:,2)),'same'))
% 
% % blind deconvolution
% [s_deconv, PSF_hat] = deconvblind(s,PSF_i);
% %%
% figure(2)
% subplot(211)
% plot([conv(s_deconv,PSF_hat,'same') s])
% legend('deconv.','signal')
% title('original vs. reconstructed');
% subplot(212)
% plot(s - conv(s_deconv,PSF_hat,'same'))
% title('extracted noise');
% % subplot(313)
% % plot(cat(1,noise, s - conv(s_deconv,PSF_hat,'same')))
% %% BDN test
% 
% k = 1;
% x = sig_mat(:,k:k+1);
% snr = 0;
% % De-noise noisy signal using soft heuristic SURE thresholding 
% % and scaled noise option, on detail coefficients obtained 
% % from the decomposition of x, at level 5 by sym8 wavelet. 
% lev = 6;
% xd = wden(x,'heursure','s','one',lev,'sym8');
% 
% % Plot signals. 
% subplot(611), plot(xd(:,1)), axis tight; 
% % title('Original signal'); 
% subplot(612), plot(x), axis tight; 
% title(['Noisy signal - Signal to noise ratio = ',... 
% num2str(fix(snr))]); 
% subplot(613), plot(xd), axis tight; 
% title('De-noised signal - heuristic SURE'); 
% 
% % De-noise noisy signal using soft SURE thresholding 
% xd = wden(x,'heursure','s','one',lev,'sym8');
% 
% % Plot signal. 
% subplot(614), plot(xd), axis tight; 
% title('De-noised signal - SURE');
% 
% % De-noise noisy signal using fixed form threshold with 
% % a single level estimation of noise standard deviation. 
% xd = wden(x,'sqtwolog','s','sln',lev,'sym8');
% 
% % Plot signal. 
% subplot(615), plot(xd), axis tight; 
% title('De-noised signal - Fixed form threshold');
% 
% % De-noise noisy signal using minimax threshold with 
% % a multiple level estimation of noise standard deviation. 
% xd = wden(x,'minimaxi','s','sln',lev,'sym8');
% 
% % Plot signal. 
% subplot(616), plot(xd), axis tight; 
% title('De-noised signal - Minimax');
% 
% % If many trials are necessary, it is better to perform 
% % decomposition once and threshold it many times:
%  
% % decomposition. 
% [c,l] = wavedec(x,lev,'sym8');
%  
% % threshold the decomposition structure [c,l].
% xd = wden(c,l,'minimaxi','s','sln',lev,'sym8');
% xd = reshape(xd,size(x,1),[]);
% subplot(611), plot(xd(:,1)), axis tight; 
% subplot(612), plot(xd(:,2)), axis tight; 