function [ CP, mean_FR, dt_bin ] = CountingProcess( time, spiketime ,new_SR ) 
% =============================================== CountingProcess =====================================================================
% GENERAL : 
% Here we build an counting process. 
% CASE 2 adaptive binning. each AC gets a optimal dt. 
% CASE 3 The sampling rate will be 500 Hz
%     Input: time 
%            spiketime - spikes time
%            new_SR - new sampling rate for CPdownSampleing
%     Output: CP - CountigProcess 
%             mean_FR
%             dt_bin
% =====================================================================================================================================

sortedIntervals = sort(diff(spiketime));

    switch nargin
        case 2 % dt is optmial and specific per CP
            dt_bin = sortedIntervals(round( 0.1 * length( spiketime )/100 ) + 1); % new dt considers only 95% of the spike as seperated spikes
            time_binned = linspace(0,time(end),time(end)/(dt_bin)); % new time axis
            mean_FR = (length(spiketime)/(length(time_binned)*dt_bin)).^2; % expexted mean firing rate ^2   
                  
            % Small value of E(lambda)^2 is biologically and comput problematic ( poor impotent acorr) the next iterative loop fixes this problem
            while mean_FR < 1.1
                % the 0.857 parameter is empirically optimized for The
                % problem of an AC with E2_lam_est lower than 1.1 (not reasonable)
                dt_bin = 0.857 * dt_bin; 
                mean_FR = (length(spiketime)/(length(time_binned)*dt_bin)).^2;
            end
            
            time_binned = linspace(0,time(end),time(end)/dt_bin);
            CP = histc(spiketime,time_binned);

        case 3 % dt is constant for all CP
            CP = histc(spiketime,time);
            CP = CPdownSampling(CP , 1e4/new_SR ); % 10kHz to new_SR Hz
            dt_bin  = 1/new_SR;
            time_binned = linspace(0,time(end),time(end)/(dt_bin)); % new time axis
            mean_FR = (length(spiketime)/(length(time_binned)*dt_bin)); %mean firing rate expectation
    end
    
end