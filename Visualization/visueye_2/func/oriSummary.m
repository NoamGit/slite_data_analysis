function [log_table, tuning_curve_data] = oriSummary(handles,cell_data, param)
%oriSummary(handles,cell_data); returns polar plots and anova tests
% to input details

if(nargin < 3)
    lw = 1; % linewidth
else
    if isfield(param,'lineWidth'); lw = param.lineWidth;else lw =1; end
    if isfield(param,'baseline'); baseline = param.baseline;else baseline = 0; end
    if isfield(param,'handles'); handles = param.handles;else handles = 0; end
end

l = cell_data.props{1}.signalLenght;
if(isstruct(cell_data.stim{:}))
    cell_data.stim = {cell_data.stim{:}.stim};
end

% quantize stimulus
stim = quantizeOriStim(cell_data);

% find off response stimulus
stim_off = calcOffResponseStim(stim);

[all_data_S, all_data_C, all_data_F] = deal(zeros(1,l));
all_baseline = 0;
num_reps = numel(cell_data(:,1));

mean_s = cell_data.data.S;
mean_c = cell_data.data.C;
mean_f = cell_data.data.Df;
all_baseline = baseline;

tidx = isnan(mean_s);
mean_s(tidx) = 0;
polar_vec = accumarray( stim',mean_s,[],@mean );
mean_s(tidx) = NaN;
polar_vec = [polar_vec(7:end)' polar_vec(2:6)'];
polar_vec = polar_vec./max(polar_vec);

% polar plot
theta = linspace(0, 2*pi, 9);
rho = [polar_vec';polar_vec(1)];

% Polar plots
fig2 = figure(2);clf;
ax2 = axes;
p = polar(theta,rho');
all_lines = findall(gca,'type','line'); 
polarticks(8,all_lines(1), ax2);
set(ax2,'NextPlot','replacechildren');
set(fig2,'position',[ 300 500 450 410]);

% set(ax2,'Box','off','visible','off');
cla(handles.axes5)
set(handles.axes5,'Visible','off');
axes(handles.axes5);
set(handles.axes5,'NextPlot','replacechildren');
new_handle = copyobj(allchild(ax2),handles.axes5);

% muliple comparison 1-way Anova test
% TODO:: make subpart/substate of 3 seconds after stimulus
% disappears.Otherwise Anova won't indicate on significance. try these line
% of codes
% deg_vec(stim ==1 & mean_s == 0) = [];
% mean_s(stim ==1 & mean_s == 0) = [];

mean_s(tidx) = eps;
% degree_map = {'blank','135 deg','180 deg','225 deg','270 deg','315 deg','0 deg','45 deg','90 deg'};
degree_map = {'blank','135 deg','180 deg','225 deg','270 deg','315 deg','0 deg','45 deg','90 deg'};
deg_vec = degree_map(stim);
[p,t,stats] = anova1(mean_s,deg_vec,'off');
fig = figure(1);
[c,m,h,nms] = multcompare(stats,'Alpha',0.05);
tuning_curve_data.data = m;
tuning_curve_data.label = nms;
ori_1 = nms(c(:,1));
ori_2 = nms(c(:,2));
p_val = c(:,6);
log_table = table(ori_1, ori_2,p_val);
set(h,'NextPlot','replacechildren');
axes(handles.axes2);
set(gca,'NextPlot','replacechildren');
cla(handles.axes2);
copyobj(allchild(get(h,'children')),handles.axes2);
title('Multiple comparison of mean p = 0.05');
axes(handles.axes2);
xlabel('Mean of response to orientation');
set(handles.axes2,'yTicklabel',flipud([{''};nms;{''}]));
set(fig,'position',[ 800 500 450 410]);

%% for offset orentation bar
polar_vec_off = accumarray( stim_off',mean_s,[],@mean );
polar_vec_off = [polar_vec_off(7:end)' polar_vec_off(2:6)'];
polar_vec_off = polar_vec_off./max(polar_vec_off);

% polar plot
rho_off = [polar_vec_off';polar_vec_off(1)];

% Polar plots
fig10 = figure(10);clf;
ax10 = axes;
p_off = polar(theta,rho_off');
all_lines = findall(gca,'type','line'); 
polarticks(8,all_lines(1), ax10);
set(fig10,'position',[ 300 50 450 410]);

% anova
degree_map = {'blank','135 deg','180 deg','225 deg','270 deg','315 deg','0 deg','45 deg','90 deg'};
deg_vec_off = degree_map(stim_off);
[p_off,t_off,stats_off] = anova1(mean_s,deg_vec_off,'off');
fig11 = figure(11);ax11 = axes;
[c_off,m_off,h_off,nms_off] = multcompare(stats_off,'Alpha',0.05);
off_tuning_curve_data.data = m_off;
off_tuning_curve_data.label = nms_off;
ori_1 = nms_off(c(:,1));
ori_2 = nms_off(c(:,2));
p_val = c(:,6);
log_table_off = table(ori_1, ori_2,p_val);
title('Multiple off response comparison of mean p = 0.05');
xlabel('Mean of response to orientation');
nms_off{1} = 'onset';
set(gca,'yTicklabel',flipud(nms_off));
set(fig11,'position',[ 800 50 450 410]);

% OFF tuning curve
fig12 = figure(12);clf;
[labels, imap] = sortDegrees(off_tuning_curve_data.label);
[~,where_blank] = ismember('blank',off_tuning_curve_data.label);
[~,where_0] = ismember('0 deg',off_tuning_curve_data.label);
imap(imap==where_blank) = [];
x = off_tuning_curve_data.data(imap,1);
x = [x; x(1)];labels(end) = {'360 deg'};
imap =[imap;where_0];
degs = (0:pi/4:pi*2)';
[~,imax] = max(mean(reshape([x; x; x],3,[]),1));
tuning_cntr = mod(imax,numel(x));
if(tuning_cntr == 0);tuning_cntr = 8; end;

[gof,gof1,gof2] = deal(struct('sse',0,'rsquare',0,'dfe',0,'adjrsquare',0,'rmse',0));
shift = 0;
gauss_index = 1; 
[f1,f2] = deal(0);
gauss_fit = {'gauss1','gauss2'};
rand_tuning = 1;
if(~rand_tuning)
    disp('finding best shift for tuning curve...');
    for k = -5:5 % fixme
        pusher_itr = [3 0 -3 2 -1 -4 1 -2]+k; % TODO: try some offsets
        bm_itr = pusher_itr(tuning_cntr);% bookmark
        x_itr = circshift(x,[bm_itr 0]);
        try
            [f1,gof1] = fit(degs,x_itr,'gauss1');
        catch
        end
        try
            [f2,gof2] = fit(degs,x_itr,'gauss2');
        catch
        end
        all_f = {f1,f2};
        all_gof = [gof1, gof2];
        [gof_itr,gauss_index_itr] = max([gof1.rsquare, gof2.rsquare]);
        if(gof.rsquare < gof_itr)
            gauss_index = gauss_index_itr;
            gof = all_gof(gauss_index);
            shift = k;
            f = all_f{gauss_index};
        end
    end
    else
      gauss_index = 1;
      shift = 0;
end
 %
pusher = [3 0 -3 2 -1 -4 1 -2]+shift; % TODO: try some offsets
bm = pusher(tuning_cntr);% bookmark
imap_optimal = circshift(imap,[bm 0]);
x_optimal = circshift(x,[bm 0]);
[f,gof] = fit(degs,x_optimal,gauss_fit{gauss_index});

degs_interp = linspace(degs(1),degs(end),200);
cfit = feval(f,degs_interp);

% p = polyfit(degs,x_optimal,1);
% cfit = polyval(p,degs_interp);

std_bar = off_tuning_curve_data.data(imap_optimal,2);
if(min(x_optimal-std_bar) < 0 )
    [~,ind] = min(x_optimal-std_bar);
    htun = errorbar(degs,x_optimal-x_optimal(ind)+std_bar(ind),std_bar,' o');
    hold on;
    plot(degs_interp, cfit-x_optimal(ind)+std_bar(ind));
    hold off;
else
    htun = errorbar(degs,x_optimal,std_bar,' o');
    hold on;
    plot(degs_interp, cfit);
    hold off;
end

set(gca,'XTick',degs,'XTickLabel',[circshift(labels,[bm 0]);{''}],'XTickLabelRotation',90);            
title('OFF response tuning curve','FontSize',14)
ylabel('Spike counts','FontSize',11);
set(fig12,'position',[ 1200 50 450 410]);
% axis_x = get(gca,'XLim');
% axis([axis_x(1) axis_x(end) 0 0.7]);
end

function [y, map] = sortDegrees(x)
    x_num = zeros(numel(x),1);
    for c = 1:numel(x)
        x_scan = textscan(x{c},'%s','delimiter',' ');
        x_num(c) = str2double(cell2mat(x_scan{1}(1)));
    end
    [val,map] = sort(x_num);
    y = cellfun(@(c) strjoin([{num2str(c)},{'deg'}],' '), num2cell(val(1:end-1)),'UniformOutput',false);
    y = [y;{'blank'}];
end
% copyobj(allchild(allchild(h)),handles.axes2);
% close(fig);
% set(gca,'ActivePositionProperty','outerposition')
% set(gca,'Units','normalized')
% set(gca,'OuterPosition',[0 0 1 1])
% set(gca,'position',[0.1300 0.1100 0.7750 0.8150])