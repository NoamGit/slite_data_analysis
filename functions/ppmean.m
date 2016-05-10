function Ex=ppmean(times)

% Ex=ppmean(times) estimates the mean rates of the point processes
% Ex - a scalar or a cell array of estimated mesn rates
% times - a vector or a cell array of the event times of the processes

vFlag=false;
if ~iscell(times)
    vFlag=true;
    % reshape times vector to be a column vector 
    % (also returns an error if times is not a 1-D vector)
    times=reshape(times, length(times), 1);
    times=mat2cell(times, size(times, 1), size(times, 2));
end

nCells=length(times);
for iCell=1:nCells
    Ex{iCell}=length(times{iCell})/(max(times{iCell}));
end

% return a scalar (and not a cell array) if times was a vector
if vFlag
    Ex=Ex{1};
end