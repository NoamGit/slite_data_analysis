% feature analysis script
% data structure is [qi, mi, csd]

cd('C:\Users\noambox\Documents\Sync');
load('.\Neural data\stat_features\all_feat_spectral_last.mat');
load('.\Neural data\ori_table.mat');
global ori_features_filt
global lbl
global ori_table
global ori_features
%% annotate each data point with a label

fig = figure;
ori_features_filt = ori_features(ori_features.label == 1,:);
lbl = [ori_features_filt.datel ,ori_features_filt.filename];
l = length(ori_features_filt.QI);
c = (1+z1( ori_features_filt.QI )) .*255;
h0 = scatter3(ori_features_filt.QI,ori_features_filt.MI,ori_features_filt.CSD, 30,c);
h0.MarkerFaceColor = [0 .70 .70];

xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('ORI features')
set(gca,'FontSize',11);
dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn',@oriInteractive)
%%
% k = 3;
% cd 'C:\Users\noambox\Documents\Sync\Thesis writing\figs\chapter - results\Summary Plots\ORI'
% set(0, 'ShowHiddenHandles', 'on')
% now,
% all_fig = get(0, 'Children');

saveas(all_fig(4),[num2str(k),'_polar','.jpg'])
savefig(all_fig(1),[num2str(k),'_scatter','.fig'])
set(all_fig(2),'Position',[1 1 1200 400]);
set(all_fig(2),'WindowStyle','normal');
set(all_fig(2),'Position',[1 1 1200 400]);
savefig(all_fig(2),[num2str(k),'_signal','.fig'])
% saveas(all_fig(2),['signal_',num2str(k),'.png'])
savefig(all_fig(3),[num2str(k),'_anova','.fig'])
k = k +1;