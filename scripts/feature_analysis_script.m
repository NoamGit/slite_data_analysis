% feature analysis script
% data structure is [qi, mi, csd]
if(isdir('C:\Users\noambox\Documents\Sync'))
    cd('C:\Users\noambox\Documents\Sync');
else
    isdir('C:\Users\Noam\Documents\Sync')
end
load('.\Neural data\stat_features\all_feat_spectral_last.mat');
load('.\Neural data\ori_table.mat');
global ori_features_filt
global lbl
global ori_table
global ori_features
%% plot features script

ori_indx = find(prod(ori_features,2));
crp_indx = find(prod(crp_features,2));

[x1, y1, z1] = deal(ori_features(ori_indx,1),ori_features(ori_indx,2),ori_features(ori_indx,3));
figure();
h1 = scatter3(x1,y1,z1,30,color(1:length(x)));
h1.MarkerFaceColor = [0 .75 .75];
xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('ORI features')
set(gca,'FontSize',11);

[x2, y2, z2] = deal(crp_features(crp_indx,1),crp_features(crp_indx,2),crp_features(crp_indx,3));
figure();
h2 = scatter3(x2,y2,z2, 30,color(1:length(x2)));
% h2 = scatter3(x2,y2,z2,'MarkerFaceColor',[0 .75 .75]);
h2.MarkerFaceColor = [0 .75 .75];
xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('CRP features')
set(gca,'FontSize',11);

% since we have some QI information for non-repetative response we filter
% them out
ori_indx2 = find((ori_features(:,4) > 0.33) & prod(ori_features,4:6) > 0 );
crp_indx2 = find((crp_features(:,4) > 0.2) & prod(crp_features,4:6) > 0 );
%% for python labeling ORI
ori_features_array = table2array( ori_features(:,(4:7)));
ori_features_array((ori_features_array(:,1) < 0.33 |...
    prod(ori_features_array(:,1:3),2) == 0 ),:) = [];
indx = (1:size(ori_features_array,1));
h1 = scatter3(ori_features_array(:,1),ori_features_array(:,2),...
    ori_features_array(:,3), 15, ori_features_array(:,4));
% h2 = scatter3(x2,y2,z2,'MarkerFaceColor',[0 .75 .75]);
h1.MarkerFaceColor = [0 .90 .90];
cdata = get(h1, 'CData');
xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('Feature space for responses to ORI')
set(gca,'FontSize',14);

%% using k-means to find infor-noninfo clusters
ori_indx2 = find((ori_features(:,4) > 0.33) & prod(ori_features,4:6) > 0 );

[x1, y1, z1] = deal(ori_features(ori_indx2,1) ,ori_features(ori_indx2,2),ori_features(ori_indx2,3));
% [x1, y1, z1] = deal(norm_nc( ori_features(ori_indx2,1) ,6) ,norm_nc(ori_features(ori_indx2,2), 6),norm_nc(ori_features(ori_indx2,3),6));
c = kmeans(ori_features(ori_indx2,:,:)',2);
figure();
h3 = scatter3(x1,y1,z1, 30,c);
% h2 = scatter3(x2,y2,z2,'MarkerFaceColor',[0 .75 .75]);
h3.MarkerFaceColor = [0 .85 .85];
xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('ORI features')
set(gca,'FontSize',11);

[x2, y2, z2] = deal(crp_features(crp_indx2,1),crp_features(crp_indx2,2),crp_features(crp_indx2,3));
c = kmeans(crp_features(crp_indx2,:,:)',2);
figure();
h3 = scatter3(x2,y2,z2, 30,c);
% h2 = scatter3(x2,y2,z2,'MarkerFaceColor',[0 .75 .75]);
h3.MarkerFaceColor = [0 .85 .85];
xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('CRP features')
set(gca,'FontSize',11);

%%

h2 = scatter3(ori_features(~indx,1),ori_features(~indx,2),ori_features(~indx,3), 30,ori_features(~indx,4));
% h2 = scatter3(x2,y2,z2,'MarkerFaceColor',[0 .75 .75]);
h2.MarkerFaceColor = [0 .70 .70];
xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('CRP features')
set(gca,'FontSize',11);
 %% for python labeling CRP
 
% indx = (crp_features(:,4) == 0) | (crp_features(:,4) == 3);
indx = (1:length(crp_features));
h3 = scatter3(crp_features(indx,1),crp_features(indx,2),crp_features(indx,3), 30,crp_features(indx,4));
% h2 = scatter3(x2,y2,z2,'MarkerFaceColor',[0 .75 .75]);
h3.MarkerFaceColor = [0 .70 .70];
xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('CRP features')
set(gca,'FontSize',11);

h4 = scatter3(crp_features(~indx,1),crp_features(~indx,2),crp_features(~indx,3), 30,crp_features(~indx,4));
% h2 = scatter3(x2,y2,z2,'MarkerFaceColor',[0 .75 .75]);
h4.MarkerFaceColor = [0 .70 .70];
xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('CRP features')
set(gca,'FontSize',11);

%% summery plots

figure()
h1 = scatter3(temp.QI,temp.MI,temp.CSD, 30);
% h2 = scatter3(x2,y2,z2,'MarkerFaceColor',[0 .75 .75]);
h1.MarkerFaceColor = [0 .70 .70];
xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('ORI features')
set(gca,'FontSize',11);

figure()
h2 = scatter3(crp_features(:,1),crp_features(:,2),crp_features(:,3), 30,crp_features(:,4));
% h2 = scatter3(x2,y2,z2,'MarkerFaceColor',[0 .75 .75]);
h2.MarkerFaceColor = [0 .70 .70];
xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('CRP features')
set(gca,'FontSize',11);
%% last try

ori_indx2 = find((table2array(ori_features(:,4)) > 0.33) & prod(table2array(ori_features(:,4:6)),2) > 0 );
crp_indx2 = find((table2array(crp_features(:,4)) > 0.2) & prod(table2array(crp_features(:,4:6)),2) > 0 );
%% annotate each data point with a label

fig = figure;
ori_features_filt = ori_features(ori_features.label == 1,:);
% lbl = strcat(temp.datel,strcat({'_'},temp.filename));
% lbl = temp.datel;
lbl = [ori_features_filt.datel ,ori_features_filt.filename];
l = length(ori_features_filt.QI);
c = (1+z1( ori_features_filt.QI )) .*255;
h0 = scatter3(ori_features_filt.QI,ori_features_filt.MI,ori_features_filt.CSD, 30,c);
h0.MarkerFaceColor = [0 .70 .70];
% a = [1:l]'; b = num2str(a); c = cellstr(b);
% dx = 0.001; dy = 0.001;dz = 0.001; % displacement so the text does not overlay the data points
% text(temp.QI+dx,temp.MI+dy,temp.CSD+dz, lbl);

xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('ORI features')
set(gca,'FontSize',11);
dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn',@oriInteractive)
%% figure for thesis section results - ORI Data filtering

figure;
axes1 = axes('Parent',gcf,'FontSize',11);
view(axes1,[-16.5 48]);
grid(axes1,'on');
hold(axes1,'on');
% Create scatter3
r_info = [ori_features.label == 1];
c(ori_features.label == 1) = 50;
scatter3(ori_features(~r_info,:).QI,ori_features(~r_info,:).MI,...
    ori_features(~r_info,:).CSD,20,'b','DisplayName','r non-informative',...
    'MarkerFaceColor',[0.39215686917305 0.474509805440903 0.635294139385223],...
    'MarkerEdgeColor',[0.152941182255745 0.227450981736183 0.372549027204514]);
% Create scatter3
scatter3(ori_features(r_info,:).QI,ori_features(r_info,:).MI,...
    ori_features(r_info,:).CSD,20,'r','DisplayName','r informative',...
    'MarkerFaceColor',[0.917647063732147 0.505882382392883 0.505882382392883],...
    'MarkerEdgeColor',[0.717647075653076 0.172549024224281 0.172549024224281]);

xlabel('QI','FontSize',24);
ylabel('MI','FontSize',24);
zlabel('CSD','FontSize',24);
title('ORI features','FontSize',24);
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.153884948931425 0.727179028934314 0.213359915544599 0.0790513812788862],...
    'FontSize',16);

disp(['number of informative r ',num2str(nnz(r_info)),' num of non-info ',num2str(nnz(~r_info))])
%% figure for thesis section results - CRP Data filtering

figure;
axes1 = axes('Parent',gcf,'FontSize',11);
view(axes1,[-16.5 48]);
grid(axes1,'on');
hold(axes1,'on');
% Create scatter3
r_info = [crp_features.label == 1];
c(crp_features.label == 1) = 50;
scatter3(crp_features(~r_info,:).QI,crp_features(~r_info,:).MI,...
    crp_features(~r_info,:).CSD,20,'b','DisplayName','r non-informative',...
    'MarkerFaceColor',[0.39215686917305 0.474509805440903 0.635294139385223],...
    'MarkerEdgeColor',[0.152941182255745 0.227450981736183 0.372549027204514]);
% Create scatter3
scatter3(crp_features(r_info,:).QI,crp_features(r_info,:).MI,...
    crp_features(r_info,:).CSD,20,'r','DisplayName','r informative',...
    'MarkerFaceColor',[0.917647063732147 0.505882382392883 0.505882382392883],...
    'MarkerEdgeColor',[0.717647075653076 0.172549024224281 0.172549024224281]);

xlabel('QI','FontSize',24);
ylabel('MI','FontSize',24);
zlabel('CSD','FontSize',24);
title('CRP features','FontSize',24);
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.153884948931425 0.727179028934314 0.213359915544599 0.0790513812788862],...
    'FontSize',16);

disp(['number of informative r ',num2str(nnz(r_info)),' num of non-info ',num2str(nnz(~r_info))])