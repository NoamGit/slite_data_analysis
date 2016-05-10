function [ sta, w_sta ] = compute_white_sta( stimulus, CP, numSamples )
%computes STA according to LS solution
%  Description:     
%           1. finds X - the 'Spike Triggered Ensemble' ('STE' - a collection of all
%           stimulus paths that happen before a spike). 
%           2. calculates LS by k = inv(X' * X) * X * Y
%  ref - slides 6 in MCM
    
    % 1. find STE
    rawStimuli = repmat(stimulus,1,numSamples);
    mid = ceil( size(rawStimuli,1)/2 );
    shiftMatrix = fliplr( [zeros( size(rawStimuli,1) - mid,numSamples) ; eye( mid,numSamples  )] );
    shiftfun = @(A,B) conv(A,B,'same');
    rawStimuli = cell2mat( arrayfun(@(k) shiftfun(rawStimuli(:,k),shiftMatrix(:,k)),...
        (1:numSamples),'UniformOutput',false));
    X_squeeze = rawStimuli( logical(CP),:);
    Y = CP( logical(CP) );
    X = bsxfun(@times, X_squeeze ,Y);
    
    % 2. find sta 
    sta = mean(X,1)';
    w_sta = inv(X' * X) * X' * Y;
end

