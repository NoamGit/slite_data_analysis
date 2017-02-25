%% this is the final analysis scrip for the chirp neural responses

%% create data frame
global celltype
global stim_map

% load cell dictionary (according to ORI)
parent_path = pwd;
parent_path = textscan(parent_path,'%s','delimiter','\');
cpath = strjoin(parent_path{:}(1:3) ,'\');
load([cpath,'\Documents\Sync\neural data\stim_features\ori_feature_mat_thresh1.mat']);
id = table2save.Properties.RowNames;
celltype = table(table2save.onoff,table2save.on_off,'RowNames',id,...
    'VariableNames',{'ONorOFF' 'ONOFF'});
load_struct = load([cpath,'\Documents\Sync\Neural data\table_data\thres1\maps.mat']);
stim_map = load_struct.stim_map;

% open visueye_2 and load chirp analysis fun (This was done manually)
% visueye_2;

% load result df
load([cpath,'\Documents\Sync\code\Data Analysis\# chirp analysis\data\crp_last_analysis.mat']);
%% data frame analysis

% fix different length problem
fix_crp = @(x) x(:,1:60);
fix_amp = @(x) x(:,1:132);
table2save.amp_mat_s = cellfun(fix_amp, table2save.amp_mat_s,'UniformOutput',false);
table2save.amp_mat_c = cellfun(fix_amp, table2save.amp_mat_c,'UniformOutput',false);
table2save.crp_mat_s = cellfun(fix_crp, table2save.crp_mat_s,'UniformOutput',false);
table2save.crp_mat_c = cellfun(fix_crp, table2save.crp_mat_c,'UniformOutput',false);

% step zero - divide to f and s gcamp types
df_f = table2save(strcmp(table2save.gcamp_type,'f'),:);
df_s = table2save(strcmp(table2save.gcamp_type,'s'),:);

% step 1.1 - divide s into on, off, on-off and else df
df_s_on = df_s(strcmp(df_s.cell_type,'on'),:);
df_s_off = df_s(strcmp(df_s.cell_type,'off'),:);
df_s_on_off = df_s(strcmp(df_s.cell_type,'on_off'),:);

% step 1.2 - divide f into on, off, on-off and else df
df_f_on = df_f(strcmp(df_f.cell_type,'on'),:);
df_f_off = df_f(strcmp(df_f.cell_type,'off'),:);
df_f_on_off = df_f(strcmp(df_f.cell_type,'on_off'),:);
%% PLOT 2 the ON OFF comparison

% plots parameters

on_off_thres = 0.05;
nBins = 15;
gcamp = 'f';
df_gcamp = struct('on', df_f_on,'off',df_f_off,'on_off',df_f_on_off);
%% plots 2.1

    on_val_on = cell2mat(cellfun(@(c) reshape(c,[],1),df_gcamp.on.on_mat_s,'UniformOutput',false)); 
    non0_on = find(on_val_on > on_off_thres );
    on_val_off = cell2mat(cellfun(@(c) reshape(c,[],1),df_gcamp.off.on_mat_s,'UniformOutput',false));
    non0_off = find(on_val_off > on_off_thres);
    on_val_oo = cell2mat(cellfun(@(c) reshape(c,[],1),df_gcamp.on_off.on_mat_s,'UniformOutput',false));
    non0_oo = find(on_val_oo > on_off_thres);

    figure(1)
    h1 = histogram(on_val_off(non0_off),nBins,'Normalization','probability');hold on;
%     histogram(on_val_oo(non0_oo),'Normalization','probability','BinWidth',h1.BinWidth);
    histogram(on_val_on(non0_on),'Normalization','probability','BinWidth',h1.BinWidth);
    hold off;
%     legend('off cell','on-off cell','on cell');
    legend('off cell','on cell');
    xlabel('spike values');
    title(['Value distribuation ON stimulus duration with GCaMP6',gcamp]); set(gca,'FontSize',20);
%% plots 2.2

    off_val_on = cell2mat(cellfun(@(c) reshape(c,[],1),df_gcamp.on.off_mat_s,'UniformOutput',false));
    non0_on = find(off_val_on > on_off_thres);
    off_val_off = cell2mat(cellfun(@(c) reshape(c,[],1),df_gcamp.off.off_mat_s,'UniformOutput',false));  
    non0_off = find(off_val_off > on_off_thres);
    off_val_oo = cell2mat(cellfun(@(c) reshape(c,[],1),df_gcamp.on_off.off_mat_s,'UniformOutput',false)); 
    non0_oo = find(off_val_oo > on_off_thres);

    figure(2)
    h2 = histogram(off_val_on(non0_on),nBins,'Normalization','probability');hold on;
%     histogram(off_val_oo(non0_oo),'Normalization','probability','BinWidth',h2.BinWidth);
    histogram(off_val_off(non0_off),'Normalization','probability','BinWidth',h2.BinWidth);
    hold off;
%     legend('on cell','on-off cell','off cell')
    legend('on cell','off cell')
    title(['Value distribuation OFF stimulus duration with GCaMP6',gcamp])
    xlabel('spike values');
    set(gca,'FontSize',20);
%% PLOT 3 the CHIRP comparison

    % plots 3.1
    load([cpath,'\Documents\Sync\code\Data Analysis\# chirp analysis\data\f_chirp.mat']);
    l = 2;
    f_axis = f_chirp(1:3:end-1);
    dt = 0.1;

    crp_val_on = cell2mat(cellfun(@(c) z1(mean(c,1)),df_gcamp.on.crp_mat_s,'UniformOutput',false));
    crp_val_off = cell2mat(cellfun(@(c) z1(mean(c,1)),df_gcamp.off.crp_mat_s,'UniformOutput',false));  
    crp_val_oo = cell2mat(cellfun(@(c) z1(mean(c,1)),df_gcamp.on_off.crp_mat_s,'UniformOutput',false)); 

    cell2plot = crp_val_oo; 
    cell2plot = reshape(cell2plot(~isnan(cell2plot)),[],size(cell2plot,2));
    figure(3)
    subplot(10,1,[1,2,3,4,5])
    h3 = imagesc( cell2plot );
    f_axis = round(f_axis.*1e2)./1e2;
    xticks = linspace(1, size(cell2plot, 2), numel(f_axis(1:5:end)));
    % set(gca, 'XTick', xticks, 'XTickLabel', f_axis(1:5:end));
    set(gca,'Xtick',get(gca,'Xtick'),'XTickLabel', str2double(get(gca,'XtickLabel')).*dt);
    ylabel('trial');xlabel('time [sec]');
    title(['Normalized mean respone to chirp with GCaMP6',gcamp])
    set(gca,'FontSize',15);

    subplot(10,1,[8,9,10])
    trail_mean = mean(cell2plot);
    f_axis_binned = binn(f_axis, l, @max);
    normalization_vect = reshape(histcounts(sort(f_axis), [0 sort(binn(f_axis, l, @max))']),[],1);
    ind2remove = find(~normalization_vect);
%     trail_mean_binned = binn(trail_mean, l, @sum);
    psth = @(x) binn(x, l, @sum)./(normalization_vect * dt); % PSTH with non linear normalization. the nrmalization factor is the bin duration
    trail_data = cell2mat(cellfun(psth, num2cell(cell2plot,2),'UniformOutput',false).');
    f_axis_binned(ind2remove,:) = [];
    normalization_vect(ind2remove,:) = [];
    trail_data(ind2remove,:) = [];
    psth = @(x) binn(x, l, @sum)./(normalization_vect * dt); % PSTH with non linear normalization. the nrmalization factor is the bin duratio
    [sorted,map] = sort(f_axis_binned,'ascend');
    [~, map2] = unique( sorted );
    sorting = map(map2);
    harea = shadedErrorBar( f_axis_binned(sorting).',trail_data(sorting,:).',...
        {@(x1) mean(x1,1), @(x2) std(x2)  }, '-b', 0);
%     harea = area(f_axis_binned(map(map2)),psth( trail_mean(map(map2))) );
%     set(harea,'FaceColor',[0.8039    0.8784    0.9686]); % Your alpha value is the 0.3
%     set(harea,'EdgeColor','none'); 
%     hold on; plot(f_axis_binned(map(map2)),psth(map(map2)));
%     hstem = stem(f_axis_binned(map(map2)),psth(map(map2)),'k'); 
%     set(hstem, 'LineWidth',0.8); 
%     hold off
    axis tight;grid on;
    ylabel('spikes per chirp freq');
    xlabel('chirp frequency [Hz]');
    title(['PSTH of response to chirp (binSize = ',num2str(l.*dt),' sec)' ])
    set(gca,'FontSize',15);
%% PLOT 4 the AMP comparison
    % plots 4.1

    amp_val_on = cell2mat(cellfun(@(c) z1(mean(c,1)),df_gcamp.on.amp_mat_s,'UniformOutput',false));
    amp_val_off = cell2mat(cellfun(@(c) z1(mean(c,1)),df_gcamp.off.amp_mat_s,'UniformOutput',false));  
    amp_val_oo = cell2mat(cellfun(@(c) z1(mean(c,1)),df_gcamp.on_off.amp_mat_s,'UniformOutput',false)); 

    cell2plot = amp_val_off; 
    cell2plot = reshape(cell2plot(~isnan(cell2plot)),[],size(cell2plot,2));
    figure(4)
    subplot(10,1,[1,2,3,4,5])
    h4 = imagesc( cell2plot );
    ylabel('trial');xlabel('time [sec]');
    title(['Normalized mean respone to AMP with GCaMP6',gcamp])
    set(gca,'FontSize',20);
    
    clear S_xx ;clear F_x;
    for k = 1:size(cell2plot,1)
        [S_xx(:,k),fs, F_x(:,k)] = ZOHSpectralEst(cell2plot(k,:),dt);
    end
    Hpsd = dspdata.psd(mean(S_xx,2),'Fs',fs);
    subplot(10,1,[8,9,10])
    plot(Hpsd);
    axis tight
    title(['PSD of response to amplitude modulation' ])
    set(gca,'FontSize',20);
    
    % plot phase
    figure(6)
    phs_on = unwrap(angle(mean(F_on,2)));
    phs_off = unwrap(angle(mean(F_off,2)));
    ly = length(mean(F_on,2));
    f = (0:ly-1)/ly*fs;

    hphase = plot(f,phs_on/pi ,f,phs_off/pi);
    set(hphase,'LineWidth' , 3)
    xlabel 'Frequency (Hz)'
    ylabel 'Phase / \pi'
    title('Phase comparison of response to amp modulation')
    legend('ON cell','OFF cell','Location','southwest')
    grid
    set(gca,'FontSize',15);
