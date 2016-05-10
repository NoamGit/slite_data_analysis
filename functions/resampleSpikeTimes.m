function [ Fs_out,spiketimes_out ] = resampleSpikeTimes( Fs, spiketimes, newrate )
%resamplePointProcess(Fs, spiketimes, stimRate ,stimulus , newrate) Since
% Neuronal activity barley cross the 100hz we want to resample the data to
% a more sparse representation. the default is 400 hz 
% important: spiketimes should come in seconds unit
% ADD: RESAMPLE OF THE STIMULUS
    
    if newrate > 500
        newrate = 500; % [hz]
        disp('The specified rate is higher than 500 hz');
    end
    
    for n = 1:numel(spiketimes)
         badIND_sr = diff(spiketimes{n}) < 2/newrate; % every activity under half of SR is supressed
         spikesSurvivalRatio = sum( ~badIND_sr )/numel(badIND_sr);
         badIND_neg = spiketimes{n} < 0; % get rid of spikes that happend before stimulus
         badIND = badIND_neg | [false; badIND_sr];
         
         if spikesSurvivalRatio < 0.95 % if less than 95% of spikes survive
             warning('PointProcessData:resamplePointProcess',['less than 95%% of spikes in process %d survived the resampling.\n'... 
                 'SR wasn''t changed and Fs_out and spiketimes_out wherent assigned.'], n);
             spiketimes_out = spiketimes;
             Fs_out = Fs;
             return;
         end
         
        spiketimes{n}( badIND ) = [];
        
    end
    spiketimes_out = spiketimes;
    Fs_out = newrate; %  Hz
end

