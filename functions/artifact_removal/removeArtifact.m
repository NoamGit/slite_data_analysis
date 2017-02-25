function [ sig_table_out, artifact_table_out ] = removeArtifact( sig_table, fName, param )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% full_field_indx
% param - 
%     lev
%     showFlag

global logvar;

% initializations
if(nargin < 3)
    lev = 10;             % wavelets levels
    showFlag     = 0;     % plots
    fs = 10;
    off_seconds = 15;
    pathName = null;
else
    if isfield(param,'wavelets_levels'); lev = param.wavelets_levels;end
    if isfield(param,'showFlag'); showFlag = param.showFlag;else showFlag = 0; end
    if isfield(param,'fs'); fs = param.fs;else fs = 10;end
    if isfield(param,'off_seconds'); off_seconds = param.off_seconds; else off_seconds = 15; end
    if isfield(param,'off_period'); off_samples = param.off_period;else off_samples = off_seconds * fs; end
    if isfield(param,'pathName'); pathName = param.pathName;else pathName = null; end
end
clear param;

full_field_s = {'CRP','WGN'};
fName_token = textscan(fName,'%s','delimiter','_');
fName_token = fName_token{1};
sig_mat = table2array( sig_table );
if(any(ismember(full_field_s,fName_token)))
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
    [dr_detrend, dr_trend] = detrend_nc( dr, param_f );
    dr_detrend = dr_detrend + abs(min(dr_detrend));

    % Remove artifact and plot results
    param.win_size = win_size;
    [ sig_mat_t, a_mat ] = ArtiRemFullField(sig_mat , dr_detrend, param);
    
    % plot
    if(showFlag)
        figure(1)
        subplot(311);plot(dr); axis tight
        subplot(312);plot(dr_trend); axis tight
        subplot(313);plot(dr_detrend); axis tight
        title('Dark region detrending');
        
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
        title('artifact removal results');
    end  
else
    %% spatio-temporal artifact removal
    
    if(any(ismember(fName_token,'BDN'))) % case where we have binary dense noise
        sig_mat_t = sig_mat; % currentely doesn't work well with the spike inference and we will set unfiltered data.
%         sig_mat_t = wden(sig_mat,'sqtwolog','s','sln',lev,'sym8');
%         sig_mat_t = reshape(sig_mat_t,size(sig_mat,1),[]);
        a_mat = [];
        
        if(showFlag)
            figure(1)
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
                title('BDN artifact removal');
            end
        end
    else if(any(ismember(fName_token ,{'CNT','ORI'})))
            
            %*********************************************%
            %**************** CONSTANTS ******************%
            %*********************************************%
            param_f.Fs_filt = 10;  % Sampling Frequency
            param_f.N     = 5;     % Order
            param_f.Fstop = 0.001;  % Stopband Frequency
            param_f.Astop = 80;    % Stopband Attenuation (dB)
            
            % load artifact
            disp('***************************************************');
            disp('stimulus is spatially changing please load artifact');
            disp('***************************************************');
            if(~pathName)
                [fName, pathName] = uigetfile('D:\# Projects (Noam)\# SLITE\# DATA\# FINAL ANALYSIS\*.xlsx','Select Cell''s Data Sheet');
            else
                stim_type = fName_token{ismember(fName_token ,{'CNT','ORI'})};
                fName_art = [fName,'Artif_'];
                fName_art = fName_art([1:9,end-5:end,10:end-6]);
                try
                    sheet_art = readtable([pathName fName_art]);
                    arti_mat = sheet_art(off_samples:end,(3:end));
                    arti_mat = arti_mat{:,:};
                    assert(size(arti_mat,1) == size(sig_mat,1));
                catch ME
                  logvar = {logvar;[pathName fName_art]};
                  sig_table_out = array2table(sig_mat,'VariableNames',sig_table.Properties.VariableNames);
                  artifact_table_out = array2table([]);
                  return;
                end
            
            % detrend art
            [arti_mat_d,arti_trend] = detrend_nc(arti_mat , param_f);
 
            if(any(ismember(fName_token,['CNT'])))
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
                if(~isempty(strfind(pathName,'160314')) || ~isempty(strfind(pathName,'160315')))
                    param.flip_stim = 1;
                end
                [ sig_mat_t, a_mat ] = ArtiRemOrientation(sig_mat, arti_mat_d, param);
            end
            
           
            % plot
            if(showFlag)
                figure(1)
                k = 4;
                subplot(311);plot(arti_mat(:,k)); axis tight
                subplot(312);plot(arti_trend(:,k)); axis tight
                subplot(313);plot(arti_mat(:,k) - arti_trend(:,k)); axis tight
                title('spatial artifact detrending');
                
                figure(2)
                offset = 1;
                for z = offset + (1:4)
                    subplot(4,1,z - offset)
                    plot([ z1(sig_mat(:,z - offset)) ]); hold on;
                    h1 = plot(z1(arti_mat_d(:,z - offset))); h1.Color(4) = 0.5;
                    hold off
                    ylabel('norm signal'); xlabel('time [sec]');
                    legend('signal', 'artifact');
                    axis tight
                end
                title('spatial signal vs. artifact');
                
                figure(3)
                range = (1:size(sig_mat_t,1)); % (2e3:3.5e3)
                for z = offset + (1:4)
                    subplot(4,1,z - offset)
                    h1 = plot(sig_mat(range,z - offset));
                    hold on; plot(sig_mat_t(range,z - offset)); h1.Color(4) = 0.5;
                    h2 = plot(-1 * a_mat(range,z - offset) + min(sig_mat_t(range,z - offset)),'k'); h2.Color(4) = 0.6;
                    ylabel('F_t'); xlabel('time [sec]');
                    legend('s','s - a',' - a');
                    hold off;
                    axis tight
                end
                title('spatial artifact removal');
            end
            
        end
    end
end

sig_table_out = array2table(sig_mat_t,'VariableNames',sig_table.Properties.VariableNames);
if(isempty(a_mat))
    artifact_table_out = array2table(a_mat);
else
    artifact_table_out = array2table(a_mat,'VariableNames',sig_table.Properties.VariableNames);
end

end

