function [ feature_struct ] = cnt_feature_extraction( cell_data, h )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
jump = 3;
alpha = 0.05;
feature_struct = [];

% get summary features for cell
[ mean_response, keys_cpd, max_contrast,cell_data, group_indx ] = cntSummary(h, cell_data); 

% ANOVA testing for significane
group_cpd = group_indx(:,2);
new_grouping = reGroup( group_cpd, jump );
[~,~,stats] = anova1(cell_data.data.S,new_grouping,'off');
[c,m,~,~] = multcompare(stats,'Alpha',0.05,'Display','off');

% if significant add tag group and add features
if any(c(:,end)<alpha)
    [~,label] = max(m(:,1));
    feature_struct = struct('mean_response',{{mean_response}}, 'cpd_keys',{keys_cpd},'label',label);
else
    label = 4;
    feature_struct = struct('mean_response',{{mean_response}}, 'cpd_keys',{keys_cpd},'label',label);
end

end

function [new_grouping] = reGroup( group, jump )
    [val,~,ic] = unique(reshape(group,[],1));
    new_val = val(1:jump:end);
    new_val = reshape(repmat(new_val,1,jump)',[],1);
    new_val = new_val(1:length(val));    
    new_grouping = new_val(ic);
end

