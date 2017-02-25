% script 150721 Consulting Shy

%% SNR

avr = load('D:\# Projects (Noam)\# SLITE\# DATA\150721Retina - SHY\SNR\avr.txt');
sig = load('D:\# Projects (Noam)\# SLITE\# DATA\150721Retina - SHY\SNR\sig.txt');

size = linspace(0,490,length(avr));
avr = avr(:,2);  
sig = sig(:,2); 

mu_sig = mean(sig(92:157));
sigma_noise = std( [sig(1:56)' sig(177:end)'] );
SNR = mu_sig/sigma_noise;

plot(size, sig,'r', size, avr,'b')
hold on

%% traces

signal = load('D:\# Projects (Noam)\# SLITE\# DATA\150721Retina - SHY\traces\1_cell3.txt');
background = load('D:\# Projects (Noam)\# SLITE\# DATA\150721Retina - SHY\traces\1_bg.txt');
% stim = load('D:\# Projects (Noam)\# SLITE\# DATA\150630Retina - ANALYSIS\Text  - SHY\figures\stimulus values.txt');
f_sig = 20; % Hz
% f_stim = 1/5020; % Hz

span = 0.0009;
% signal = signal(5:end,:); background = background(5:end,:);
[ delF,t ] = DeltaF( signal, span, background);
% stim_RESAMPLED = resample(stim(:,2), f_sig, f_stim);
% stim_RESAMPLED = stim_RESAMPLED(1:length(delF));
figure(1)
plot(t/f_sig, delF,'b');
title('\Delta F/F Cell 1'); xlabel('time [sec]'); ylabel('\Delta F/F');
%%
signal = load('D:\# Projects (Noam)\# SLITE\# DATA\150721Retina - SHY\traces\1_cell1.txt');
background = load('D:\# Projects (Noam)\# SLITE\# DATA\150721Retina - SHY\traces\1_bg.txt');
f_sig = 10; % Hz

stimOn = 20; % On
stimOff = 5e3; % off
stim  = zeros( length(signal),1);
stim(31:32:end) = 0.5;

span = 0.0001;
% stim = 
[ delF,t ] = DeltaF( signal, span, background );
figure(1)
plot(t/f, delF, t/f, stim);
title('\Delta F/F Cell 1 Stim Rate 1Hz'); xlabel('time'); ylabel('\Delta F/F'); 
legend('response','stimulus');