%% script for 10.06.15 Analysis

%% Text
    % Calculate DeltaF/F

    signal = load('D:\# Projects (Noam)\# SLITE\# DATA\150630Retina - ANALYSIS\Text  - SHY\figures\Values_1.txt');
    background = load('D:\# Projects (Noam)\# SLITE\# DATA\150630Retina - ANALYSIS\Text  - SHY\figures\background.txt');
    stim = load('D:\# Projects (Noam)\# SLITE\# DATA\150630Retina - ANALYSIS\Text  - SHY\figures\stimulus values.txt');
    f_sig = 50; % Hz
    f_stim = 60; % Hz
    
    span = 0.0009;
    % signal = signal(5:end,:); background = background(5:end,:);
    [ delF,t ] = DeltaF( signal, span, background);
    stim_RESAMPLED = resample(stim(:,2), f_sig, f_stim);
    stim_RESAMPLED = stim_RESAMPLED(1:length(delF));
    figure(1)
    plot(t/f_sig, delF,'b',t/f_sig,normax(stim_RESAMPLED),'r');
    title('\Delta F/F Cell 1'); xlabel('time [sec]'); ylabel('\Delta F/F'); 
 %% Bar

    signal = load('D:\# Projects (Noam)\# SLITE\# DATA\150630Retina - ANALYSIS\Moving Bar\figures\Values_6.txt');
    background = load('D:\# Projects (Noam)\# SLITE\# DATA\150630Retina - ANALYSIS\Moving Bar\figures\background.txt');
    f_sig = 25; % Hz

    span = 0.002;
    % signal = signal(5:end,:); background = background(5:end,:);
    [ delF,t ] = DeltaF( signal, span, background);
    figure(2)
    plot(t/f_sig, delF);
    title('\Delta F/F Cell 2'); xlabel('time [sec]'); ylabel('\Delta F/F'); 