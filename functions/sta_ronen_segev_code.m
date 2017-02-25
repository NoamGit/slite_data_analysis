function [ STA ] = RonenSTA(spikesSAMPLES, Stim, spike_Fs)
%RonenSTA STA evalutaion according to Ronen Segev's code
%   Simple STA implementation. The lenght is constant (35) and the only
%   input is the stimulus itself. 
%   One can see that the stimulus is truncated in order to avoid spiked
%   that accure before 1790 [sec]
   
    spikesInSampels = round(spikesSAMPLES/(spike_Fs/30)); % this manipulation allows to move from sp [sampl] to W [sampl] resampling from 10^4 to 30 Hz in samples
    spikesInSampels = spikesInSampels( spikesInSampels<53700 ); % this implies that we deal with spikes which occure before W[53700] ~  1790 
    STA = zeros(35,1); %only 35 sampelings are accounted in the STA
    for i=20:(length(spikesInSampels)-40);
        STA = STA + Stim(spikesInSampels(i)+(-29:5)); %in each spike accumulate 35 values 
    end;
    STA = STA/(length(spikesInSampels)-61); % division neccesery for average
end

