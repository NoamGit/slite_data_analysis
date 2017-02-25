function R=ppcorr(times, dt, maxlag, bias, correctR0)

% R=ppcorr(times, dt, maxlag, bias, correctR0) calculates the correlation
%       structure of the point processes
% times - cell array (may be a vector ) of spike times of different units vector (for a single process) or
% dt - desired temporal resolution
% maxlag - maxlag (as in xcorr) - full length if not specified
% bias - bias flag (as in xcorr, on of: {'none', 'biased', ['unbiased']})
% correctR0 - flag whether to remove a 'delta' at Rxx(tau=0)
% (1 - remove, 0 - don't remove). Default set to 0.

% check the input apply defaults
if nargin<3
    maxlag=[];
end

if nargin<4
    bias='unbiased';
end

if nargin<5
    correctR0=0;
end

% a patch for a single vector times input
vFlag=false;
if ~iscell(times)
    vFlag=true;
    % reshape times vector to be a column vector
    % (also returns an error if times is not a 1-D vector)
    times=reshape(times, length(times), 1);
    times=mat2cell(times, size(times, 1), size(times, 2));
end


% calculate the correlation structure
nCells=length(times);
R=cell(nCells);
for iCell=1:nCells
    times{iCell}=times{iCell}-min(times{iCell});
    Tmax=max(times{iCell});
    Tmin=min(times{iCell});
    counts1=histc(times{iCell}, Tmin:dt:Tmax);
    if isempty(maxlag)
        maxlag=length(counts1);
    end
    R{iCell, iCell}=xcorr(counts1, maxlag, bias)/dt^2;

    % remove delta at tau==0
    if correctR0
        zerolag=(length(R{iCell, iCell})+1)/2;
        E=ppmean(times{iCell});
        R{iCell, iCell}(zerolag)=R{iCell, iCell}(zerolag)-E/dt; % this implementation is using the moment generating function
    end
    
    for jCell=iCell+1:nCells
        Tmax=max(max(times{iCell}), max(times{jCell}));
        Tmin=min(min(times{iCell}), min(times{jCell}));
        counts1=histc(times{iCell}, Tmin:dt:Tmax);
        counts2=histc(times{jCell}, Tmin:dt:Tmax);
        if isempty(maxlag)
            maxlag=length(counts1);
        end
        R{iCell, jCell}=xcorr(counts1, counts2, maxlag, bias)/dt^2;
        %using symmetry to save computations
        R{jCell, iCell}=flipud(R{iCell, jCell});
    end
end

% return a vector (and not a cell array) if times was a vector
if vFlag
    R=R{1};
end
