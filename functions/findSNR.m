%% find SNR of single trace
% Take Raw data and find Peak DF/F0 response over the 
% standard deviation of the signal during one second period before stimulation)
% load data of DF but without Kalmann filtering

function SNR = findSNR(cellInd)
  
    %% load data
    [fName, pathName] = uigetfile('D:\# Projects (Noam)\# SLITE\# DATA\*.xlsx','Select Cell''s Data Sheet');
    sheet = xlsread([pathName fName]);
    signalMat = sheet(:,(3:end)); 
    time = sheet(:,1);
    dt = double(time(2)-time(1));
    stim = sheet(:,2);
    %% data cleansing (artifact removal) and max peak detection
    
    validInd = ~isnan(stim);
    if(any(~validInd))
        errordlg('There are nan values in your dataset...');
        return
    end
    stim_shift = [0 ;stim(1:end-1)]; 
    stim = (~(stim_shift & stim) & stim);
    
    samples = (1:size(stim,1));
    stimSamples = samples(stim);
    X = signalMat(:,cellInd);
    afterStimLenght = 4/dt; % 2 sec
    stimSamples_begin = find(stim);
    peakSamp = repmat( stimSamples_begin, 1, afterStimLenght );
    peakSamp = bsxfun(@plus, peakSamp, (1:afterStimLenght) )';
    peakSamp(peakSamp > length(X)-1) = length(X)-1;
    peaks = max( X(peakSamp) );
    %% finding std of 1 sec before stimulus
    
    beforeStimLenght = 1/dt; % 2 sec
    blSamp = repmat( stimSamples_begin, 1, beforeStimLenght );
  	blSamp = bsxfun(@minus, blSamp, (1:beforeStimLenght) )';
    blSamp(blSamp > length(X)-1) = length(X)-1;
    bl = std( X(blSamp) );
    
    SNR = mean(peaks)/mean(bl);
end