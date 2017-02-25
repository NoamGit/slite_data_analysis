function [ feature_struct ] = ori_feature_extraction( cell_data, h )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

feature_struct = [];
size_on_off_wn = 20; % 2 sec win
alpha = 0.05;
num_of_ss_for_ds = 2; % number of stat_sig decisions to be considered DS

% initialization of fields 
[feature_struct.on_off,feature_struct.onoff, feature_struct.os,...
    feature_struct.os_type, feature_struct.ds, feature_struct.ds_type] = deal(nan);

% ANOVA test ON or OFF cell
group_on_or_off = cell_data.stim{:}.stim > 0;
group_on_wn = circshift([conv( double((diff(group_on_or_off)) > 0.5),[0 ones(1,size_on_off_wn)],'same'); 0],size_on_off_wn/2);
group_off_wn = circshift([conv( double((diff(group_on_or_off)) < -0.5),[0 ones(1,size_on_off_wn)],'same'); 0],size_on_off_wn/2);
indx = group_off_wn | group_on_wn;
group_onoff = group_off_wn(indx) + 2.*group_on_wn(indx);
group_maybe_vs_nochance = group_on_or_off | group_off_wn;
cell_data.data.S(isnan(cell_data.data.S)) = 0;
[p_on_off,~, stats] = anova1(cell_data.data.S,group_on_or_off,'off');
if p_on_off < alpha
    [~,argmax] = max(stats.means);
    feature_struct.on_off = str2double(stats.gnames(argmax));
end

% ANOVA for ON_OFF cell 2 tests
p_onoff_1 = anova1(cell_data.data.S(indx),group_onoff,'off');
p_onoff_2 = anova1(cell_data.data.S,group_maybe_vs_nochance,'off');
if p_onoff_1 >= alpha && p_onoff_2 < alpha
    feature_struct.onoff = 1;
end

% check OS and DS
stim_on = quantizeOriStim(cell_data);
if feature_struct.on_off == 1 % on cell
    indx = stim_on == 1;
    stim_on(indx) = [];
    stim = stim_on;
elseif feature_struct.on_off == 0
    stim_off = calcOffResponseStim(stim_on);
    indx = stim_off == 1;
    stim_off(indx) = [];
    stim = stim_off;
elseif feature_struct.onoff == 1
    stim_off = calcOffResponseStim(stim_on);
    stim_off(stim_off == 1) = 0;
    stim_on(stim_on == 1) = 0;
    stim_on_off =  stim_off + stim_on;
    indx = stim_on_off == 0;
    stim_on_off(indx) = [];
    stim = stim_on_off;
else
    feature_struct = [];
    return;
end

% check if OS
degree_map = {'blank','135/315 deg','0/180 deg','45/225 deg','90/270 deg','135/315 deg','0/180 deg','45/225 deg','90/270 deg'};
deg_vec = degree_map(stim);
[~,~,stats] = anova1(cell_data.data.S(~indx),deg_vec,'off');
[c_os,m_os,~,nms_os] = multcompare(stats,'Alpha',0.05,'Display','off');
ss_ind = c_os(:,6) < alpha;
% except_ori = [];
if(any(ss_ind))
%     intersec_itr = orient_comp{1};
%     for k = 2:size(orient_comp,1)
%         intersec_itr = intersect(intersec_itr,orient_comp{k});
%     end
%     except_ori = intersec_itr;
%     if ~isempty(except_ori)
    subs = reshape(cell2mat(num2cell(c_os(ss_ind,1:2),2)),[],1);
    [highest_mean, argmax] = max(m_os(unique(subs),1));
%     orient_comp = num2cell(c_os(ss_ind,1:2),2);
    feature_struct.os = nms_os(argmax);
    feature_struct.os_type = m_os(argmax,1) - mean(m_os(:,1));
end

% check if DS
degree_map = {'blank','135 deg','180 deg','225 deg','270 deg','315 deg','0 deg','45 deg','90 deg'};
deg_vec = degree_map(stim);
[~,~,stats_ds] = anova1(cell_data.data.S(~indx),deg_vec,'off');
[c_ds,m_ds,~,nms_ds] = multcompare(stats_ds,'Alpha',0.05,'Display','off');
ss = c_ds(:,6) < alpha;
subs = reshape(cell2mat(num2cell(c_ds(ss,1:2),2)),[],1);
% [differ, argmax] = max(accumarray(subs, ones(size(subs))));
tupels = num2cell(c_ds(ss,1:2),2);
for cand = unique(subs)'
    tupels_itr = sum(cell2mat(cellfun(@(c) ismember(c,cand),tupels,'UniformOutput',false)),2) == 1;
    neg_orie = cell2mat(cellfun(@(c) setdiff(c,cand),tupels(tupels_itr),'UniformOutput',false));
    [ neg_ories_data ] = getOrientation( neg_orie, nms_ds );
    [ cand_orie ] = getOrientation( cand, nms_ds );
    orie_range = (0:45:45 * floor(num_of_ss_for_ds/2));
%     neg_ories_h = mod(mod(cand_orie+180,360)+[fliplr(-orie_range(2:end)) orie_range],360)';
    neg_orie_h = mod(cand_orie+180,360);
    if any(ismember(neg_ories_data,neg_orie_h))
        feature_struct.ds = nms_ds(cand);
        feature_struct.ds_type =  m_ds(cand,1) - mean(m_ds(:,1));
        break; % TODO: nore than 1 Direction cannot be found in this structure
    end
end

disp(feature_struct);

end

function [ori_out] = getOrientation( ind_in, ori_dic )
    deg_str = cellfun(@(c) textscan(c,'%s','delimiter',' '), ori_dic(ind_in));
    deg_str = [deg_str{:}]';
    ori_out = mod(cellfun(@str2num,deg_str(:,1)),360);
end



