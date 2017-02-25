function [ Table ] = plotCrp(handles, fv_s, fv_c, visu_vec_s, cell_data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    % plots
    fig5 = figure(5);clf;
%     set(gca,'Position',[0.1300    0.1100    0.7750    0.2828]);
    hb = subplot(11,11,[(1:11:4*11+1),(45:49)]);
    hrose = subplot(10,10,[(6:11:3*11+1),(46:50)]);set(hrose,'Visible','off')
    hmulti = subplot(11,11,[(78:11:78+3*11+1),(112:112+9)]);
    
    % rose plot
    dict = containers.Map({'on';'off';'gray';'chirp';'amp'},[0;72;144;216;288]);
    fig1 = figure(1);clf;
    ax_1 = axes;
    Options.AngleNorth = 0;Options.AngleEast = 90;
    Options.nSpeeds = 5;Options.nDirections = 5;Options.nFreq = 3;
    Options.Min_Radius = 0.15;Options.height = 512;Options.width = 512;
    Options.cMap = 'summer';Options.legendType = 1;Options.axes = ax_1;
    Options.scalefactor = 0.92; Options.TitleString = 'Chirp features Visual';
    [~,~,~,~,Table] = WindRose(cellfun(@(c) dict(c),visu_vec_s(:,2)),...
        cell2mat(visu_vec_s(:,1)),Options);
    featureticks(ax_1);
    set(fig1,'position',[ 2111 600 450 410]);    
    cla(handles.axes2); set(handles.axes2,'Visible','off');axes(handles.axes2);
    set(handles.axes2,'NextPlot','replacechildren');copyobj(allchild(ax_1),handles.axes2);
    copyobj(allchild(ax_1),hrose);
    set(hrose,'FontSize',5);
    close(fig1);
    
    
    % bar plot
    cla(handles.axes5)
    axes(handles.axes5);
    set(handles.axes5,'NextPlot','replacechildren');
    bar([fv_s(1:3) mean(fv_s(4:8)) mean(fv_s(9:20))]...
        ,'FaceColor',[0.7569    0.8667    0.7765],'EdgeColor',[0.7569    0.8667    0.7765],'LineWidth',1.5);
    set(gca, 'XTickLabel',{'on';'off';'gray';'chirp';'amp'}, 'XTick',1:numel({'on';'off';'gray';'chirp';'amp'}))
    copyobj(allchild(handles.axes5),hb);
    set(hb, 'XTickLabel',{'on';'off';'gray';'chirp';'amp'}, 'XTick',1:numel({'on';'off';'gray';'chirp';'amp'}),'XTickLabelRotation',90)
    title(hb,'Spike Value Histogram','FontSize',22);
    set(hb,'FontSize',24);

    % multi signal plot
    fig3 = figure(3);clf;
    data_mat = cell2mat(cell_data.multi_data.multi_signal');
    data_mean = mean(data_mat,2);
    dt = 1/cell_data.properties(:).imaging_rate;
    time = (0:dt:(size(data_mean,1)-1) * dt)';
    data_stim = cell_data.multi_data.multi_stimulus{1};
    stim_scaled = z1(data_stim) .* 0.3 * (max(data_mat(:))+eps) + 1 * (max(data_mat(:))+eps);
    ax_multi = axes;
    hdata = plot(time,data_mat,'Color',[0.75    0.75    0.75]);
    axis_ax = get(ax_multi,'YLim');
    hold on;
    [~,plocs] = findpeaks(diff(diff(data_stim)),'MinPeakDistance',5,'MinPeakHeight',0.1);
    plocs = plocs + 1;
    hstim = plot(time,stim_scaled,'Color',[ 0.4706    0.3059    0.4471],'LineWidth',1.5);
    arrayfun(@(x) line([time(x) time(x)],[0 5],'LineStyle', '--'),plocs);
    hmean = plot(time,data_mean,'k','LineWidth',1.5);
    uistack(hdata,'top');
    uistack(hmean,'top');
    hold off;
    copyobj(allchild(ax_multi),hmulti);
    xlabel(hmulti,'Time [sec]','Interpreter','latex','FontSize',20);
    ylabel(hmulti,'$\hat{Ca^{+2}}$ signal','Interpreter','latex','FontSize',20);
    set(hmulti,'FontSize',20);
    title(hmulti,['                       OFF      ON     OFF Gray',...
        '         CHIRP         Gray'...
        '                                Amp  '...
        '                                 Gray  OFF                    ']...
        ,'Interpreter','tex','FontSize',12)
    axis(hmulti,[time(1) time(end) axis_ax(1) max(stim_scaled+eps)]);
    close(fig3);
end

