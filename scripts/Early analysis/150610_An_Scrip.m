%% script for 10.06.15 Analysis


%% Flashes 
% Calculate DeltaF/F

    signal = load('C:\Users\noambox\Desktop\Test Images - ImageJ\trace sample\Cell2.txt');
    background = min(signal(:,2)) * ones(length(signal),2);
    f = 31.25; % Hz

    stimOn = 10; % On
    stimOff = 1e3; % off
    stim  = zeros( length(signal),1);
    stim(31:32:end) = 0.5;

    span = 0.0001;
    % stim = 
    [ delF,t ] = DeltaF( signal, span, background );
    figure(1)
    plot(t/f, delF, t/f, stim);
    title('\Delta F/F Cell 1 Stim Rate 1Hz'); xlabel('time'); ylabel('\Delta F/F'); 
    legend('response','stimulus');
%% Moving Bar
    %% Sample 1

    signal = load('C:\Users\noambox\Videos\Retina150610\Slite_Bar\Values_Cell4.txt');
    background = load('C:\Users\noambox\Videos\Retina150610\Slite_Bar\Values_F0_2.txt');
    f = 20; % Hz

    span = 0.02;
    % signal = signal(5:end,:); background = background(5:end,:);
    [ delF,t ] = DeltaF( signal, span, background,  'bg_method', 'poly2' );
    figure(1)
    plot(t/f, delF);
    title('\Delta F/F Cell 1'); xlabel('time'); ylabel('\Delta F/F'); 
    %% Sample 2

    signal = load('C:\Users\noambox\Videos\Retina150610\Slite_Bar\Values_Cell3.txt');
    background = load('C:\Users\noambox\Videos\Retina150610\Slite_Bar\Values_F0_2.txt');
    f = 20; % Hz

    span = 0.02;
    % signal = signal(5:end,:); background = background(5:end,:);
    [ delF,t ] = DeltaF( signal, span, background, 'bg_method', 'poly2' );
    figure(2)
    plot(t/f, delF);
    title('\Delta F/F Cell 2'); xlabel('time'); ylabel('\Delta F/F'); 