function doc_datacursormode
% Plots graph and sets up a custom data tip update function
fig = figure;
temp = ori_features(ori_features.label == 1,:);
% lbl = strcat(temp.datel,strcat({'_'},temp.filename));
lbl = temp.datel;
l = length(temp.QI);
h1 = scatter3(temp.QI,temp.MI,temp.CSD, 30);
h1.MarkerFaceColor = [0 .70 .70];
% a = [1:l]'; b = num2str(a); c = cellstr(b);
% dx = 0.001; dy = 0.001;dz = 0.001; % displacement so the text does not overlay the data points
% text(temp.QI+dx,temp.MI+dy,temp.CSD+dz, lbl);

xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('ORI features')
set(gca,'FontSize',11);
dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn',@myupdatefcn)