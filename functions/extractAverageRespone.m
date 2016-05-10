%% Extract Average respone 

function extractAverageRespone()
    %% parameters 
    % TODO : Insert to gui
    windowTime = 3.5; 
    N    = 50;       % Order
    Fc   = 1;      % [Hz] Cutoff Frequency
    showSignalFlag = 1;
    num2disp = 6;
    
    %% load data
    [fName, pathName] = uigetfile('D:\# Projects (Noam)\# SLITE\# DATA\*.xlsx','Select Cell''s Data Sheet');
% %         Debug
%         pathName ='D:\# Projects (Noam)\# SLITE\# DATA\150720Retina - ANALYSIS\FLASH_20msON_20Hz_SLITE_1\';
%         fName = 'Data.xlsx';
    sheet = xlsread([pathName fName]);
    time = sheet(:,1);
    dt = double(time(2)-time(1));
    stim = sheet(:,2);
    signalMat = sheet(:,(3:end)); 
    numSignals = size(signalMat,2);
    
    %% data cleansing
    % delete conseq accurenses of stimulus
    validInd = ~isnan(stim);
    if(any(~validInd))
        errordlg('There are nan values in your dataset...');
        return
    end
    stim_shift = [0 ;stim(1:end-1)]; 
    stim = (~(stim_shift & stim) & stim);
    samples = (1:size(stim,1));
    stimSamples = samples(stim);
    stim = 5.* stim;
    
    % ** optional - LPF the signalMatrix 
    Fs = 1/dt;  % Sampling Frequency
    flag = 'scale';  % Sampling Flag
    % Create the window vector for the design algorithm.
    win = hamming(N+1);
    % Calculate the coefficients using the FIR1 function.
    b  = fir1(N, Fc/(Fs/2), 'low', win, flag);
    Hd = dfilt.dffir(b);
    filtMat = deal( zeros(size(signalMat)) );
    cleaningMat = zeros(size(signalMat(:,1),1),1);
    
    % computation
    for k=1:numSignals
        filtMat(:,k) = filtfilt( b, 1 , signalMat(:,k) );    
    end
    %% Average respone
    
    windowSize = round(windowTime * Fs); 
    if windowTime >= round((stimSamples(2) - stimSamples(1)) * dt)
        errordlg('Stimulus interval is shorter than the windowSize...');
        return;
    end
    regularWinSizeSampIndex = stimSamples+windowSize > size(signalMat(:,1),1);
    indx = repmat(stimSamples(~regularWinSizeSampIndex)',1,windowSize);
    indx = bsxfun(@plus, indx, (0:windowSize-1));
    cleaningMat( indx ) = true;
    dataStruct1 = filtMat( logical(repmat( cleaningMat,1,numSignals )) );
    dataStruct2 = reshape( dataStruct1,windowSize, [], numSignals);
    dataStruct_mean = mean(dataStruct2, 2);
    dataStruct_std = std(dataStruct2, 0 ,2);
    respMat_std = squeeze(dataStruct_std);
    respMat_mean = squeeze(dataStruct_mean);
    %% plots
    windowTimeArray = (0:dt:dt *(windowSize-1));
    plotAvrResponse( windowTimeArray, respMat_mean, respMat_std, num2disp );
    
    if(showSignalFlag)
        plotSignals(time, signalMat, stim, num2disp);
    end
end

%% plot functions

function [ ] = plotAvrResponse(time, mu, sigma, num2disp)

    % default values
    title_arg = 'Cell ';
    numOfSig = size(mu,2);
    cellIndx = (1:num2disp * ceil(numOfSig/num2disp));
    xlabelUser = 'time';
    ylabelUser = '\DeltaF/F';
    ylimit_bottom = min(mu(:))-min(sigma(:));
    axisUser = [time(1) time(end) ylimit_bottom 0.5];

    [fact] = factor(num2disp);
    x = prod(fact(1:ceil(length(fact)/2))); y =  prod(fact(ceil(length(fact)/2)+1:end));
    numFigures = ceil(numel( cellIndx )/num2disp);   
    for k = 1:numFigures
        figure(k);
        for n = 1:num2disp
            try
                dataNum = cellIndx(n+(k-1)*num2disp);
                s(n) = subplot(y,x,n);
                hE = errorbar(time, mu(:,dataNum), sigma(:,dataNum));
%                 set(hE, 'LineStyle', '--', 'Marker', '.', 'Color', [.1 .1 .1]);
                set(hE, 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
                    'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor' , [.7 .7 .7]);
                title(s(n),[title_arg,' ',num2str(dataNum)]);            
            catch err
%                 disp(err);
                break;
            end
        maxAxis = max(mu(:,dataNum)) + max(sigma(:,dataNum));
        minAxis = min(mu(:,dataNum)) - min(sigma(:,dataNum));
        if(maxAxis > 0.5)
            axisUser = [time(1) time(end) minAxis maxAxis];
        end
        axis( axisUser );
        xlabel( xlabelUser );
        ylabel( ylabelUser );
        end
    end
end

function [ ] = plotSignals(time, sigMat, stimulus, num2disp)

    % default values
    title_arg = 'Cell ';
    numOfSig = size(sigMat,2);
    cellIndx = (1:num2disp * ceil(numOfSig/num2disp));
    xlabelUser = 'time [sec]';
    ylabelUser = '\DeltaF/F';
    axisUser = [time(1) time(end) -inf inf];

    [fact] = factor(num2disp);
    x = prod(fact(1:ceil(length(fact)/2))); y =  prod(fact(ceil(length(fact)/2)+1:end));
    numFigures = ceil(numel( cellIndx )/num2disp);   
    for k = 1:numFigures
        figure();
        for n = 1:num2disp
            try
                dataNum = cellIndx(n+(k-1)*num2disp);
                s(n) = subplot(y,x,n);
                plot(time, stimulus, '--r', time, sigMat(:,dataNum));
%                 set(hE, 'LineStyle', '--', 'Marker', '.', 'Color', [.1 .1 .1]);
%                 set(hS, 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
%                     'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor' , [.7 .7 .7]);
                title(s(n),[title_arg,' ',num2str(dataNum)]);            
            catch err
%                 disp(err);
                break;
            end
        maxAxis = max(sigMat(:,dataNum));
        minAxis = min(sigMat(:,dataNum));
        if(maxAxis > 0.5)
            axisUser = [time(1) time(end) minAxis maxAxis];
        end
        axis( axisUser );
        xlabel( xlabelUser );
        ylabel( ylabelUser );
        end
    end
end