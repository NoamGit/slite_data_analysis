function h=plotcorr(R, dt, axisOnOff, scaling)

% h=plotcorr(R, dt) plots the multicorrelation structure defined in R
% h - handle to the figure
% R - cell array with correlation structure
% dt - time resolution of R
% axisOnOf - flag for turning the axis on/off. One of {['on'], 'off'}
% scaling - flag for scaling the y-axis.
%   One of {['auto'], 'same', 'cross', 'smart'}.
%     'auto' - automatic MATLAB scaling
%     'same' - all the axes are scaled to the same range
%     'cross' - the cross-correlations are 'same-scaled',
%               the auto-correlations are auto-scaled
%     'smart' - similar to 'same', but whithout taking into account
%               'deltas' in the auto-correltion functions

h=gcf; % handle to the current figure, creates a new figure if needed
set(gcf, 'Color', [0.95 0.95 0.95]);

if nargin<3
    axisOnOff='on';
end

if nargin<4
    scaling='auto';
end

nCells=length(R);
maxlag=(length(R{1,1})-1)/2;
zerolag=maxlag+1;
tau=dt*(-maxlag:maxlag);

%% calculating the y-axis margins for scaling
yMargins=[];
for iCell=1:nCells
    if isequal(scaling, 'same')
        yMargins=[yMargins, minmax(R{iCell, iCell}')];
    elseif isequal(scaling, 'smart')
        yMargins=[yMargins, minmax(R{iCell, iCell}([1:zerolag-1, zerolag+1:end])')];
    end
    for jCell=iCell+1:nCells
        yMargins=[yMargins, minmax(R{iCell, jCell}')];
    end
end
yMargins=minmax(yMargins);

%% plotting the correlation structure
for iCell=1:nCells
    for jCell=1:nCells
        %         subplot(nCells, nCells, (iCell-1)*nCells+jCell);
        subplot('Position', [(jCell-1)/nCells, (nCells-iCell)/nCells, 1/nCells, 1/nCells]);
        set(gca, 'OuterPosition',[(jCell-1)/nCells, (nCells-iCell)/nCells, 1/nCells, 1/nCells]);
        if iCell==jCell && isequal(scaling, 'smart')
            plot(tau(1:zerolag-1), R{iCell, jCell}(1:zerolag-1), '-');
            hold on;
            plot(tau(zerolag-1:zerolag+1), R{iCell, jCell}(zerolag-1:zerolag+1), ':')
            plot(tau(zerolag+1:end), R{iCell, jCell}(zerolag+1:end), '-');
        else
            plot(tau, R{iCell, jCell});
        end
        hold on;
        xlim(minmax(tau));
        if iCell==jCell
            if ismember(scaling, {'same', 'smart'})
                ylim([min(yMargins(:,1)), max(yMargins(:,2))].*[0.95, 1.05]);
            end
        else
            if ~isequal(scaling, 'auto')
                ylim([min(yMargins(:,1)), max(yMargins(:,2))].*[0.95, 1.05]);
            end
        end
        set(gca, 'Visible', axisOnOff);
        xlabel('\tau [sec]');
        ylabel(['R_{', num2str(iCell), num2str(jCell),'}']);
    end
end
