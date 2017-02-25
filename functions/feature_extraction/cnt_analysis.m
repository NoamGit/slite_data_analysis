results = table2save;
group_label = 4;
foo = cellfun(@(C)  cell2mat( cellfun(@(c) norm_nc(c,5),num2cell(C,1),'UniformOutput',false) ) ,...
        results(results.label == group_label,1).mean_response, 'UniformOutput', false);
foo = reshape(cell2mat(foo)',7,6,[]);
foo(isnan(foo)) = 0;
mean_response =  mean(foo,3)';
std_response = std(foo,[],3)';
disp(['num of cells - ',num2str(nnz(results.label == group_label))])
%% try all

keys_cpd = results.cpd_keys(1,:);
out = plotCnt([], mean_response, keys_cpd, std_response );
%% try one by one

group_label = 4;
tbl = results(results.label == group_label,:);
mean_r_vec = tbl.mean_response;
l = numel(mean_r_vec);
for k = 1:l
    out_itr{group_label}{k} = plotCnt([], mean_r_vec{k}, keys_cpd );
    contrast_thresh{group_label}{k} = out_itr{group_label}{k}.thres_contrast;
end
%%
[thres_contrast,thres_std] = deal([]);
for m = 1:numel(contrast_thresh)
    cnt_thrs = cell2mat(contrast_thresh{m});
    for k = 1:size(cnt_thrs,1)
        thres_contrast(m,k) = mean(cnt_thrs(k,~isnan(cnt_thrs(k,:))));
        thres_std(m,k) = std(cnt_thrs(k,~isnan(cnt_thrs(k,:))));
    end
% figure(6);errorbar(keys_cpd(valid_indx),thres_contrast(valid_indx),thres_std(valid_indx))
end

% averaging the responses of high-cpd responseve cells (group 3 and 4)
thres_contrast(2,:) = mean(thres_contrast(2:3,:),1);
thres_contrast(3,:) = [];

thres_contrast_mu = reshape(mean(thres_contrast,1),[],1);
valid_indx = ~isnan(thres_contrast_mu);
%% final plot

figure(2);clf;cnt_ax = axes();
% spline_thres_contrast = spline(keys_cpd(valid_indx),10.*1./thres_contrast(valid_indx),(min(keys_cpd):0.01:max(keys_cpd)));
poly_cnt = polyfit(keys_cpd(valid_indx)',10.*1./thres_contrast_mu(valid_indx),3);
poly_thres_contrast = polyval(poly_cnt,(keys_cpd(1):0.01:keys_cpd(end)));
h_cnt = loglog(cnt_ax,keys_cpd,10.*1./thres_contrast_mu,'o');
h_cnt.MarkerSize = 10;
h_cnt.MarkerFaceColor = [0.75,0.75,0.75];
grid on
%     axis([0 0.7 10 50])
% hold on;plot(cnt_ax,(keys_cpd(1):0.01:keys_cpd(end)), spline_thres_contrast); hold off
hold on;
h_plt0 = plot(cnt_ax,(keys_cpd(1):0.01:keys_cpd(end)), poly_thres_contrast); 
set(h_plt0,'LineWidth',4);
hold off
title('Contrast sensitivity function')
xlabel('Spatial freq. [cpd]');
ylabel('Sensetivity');
set(gca,'FontSize',20);

hold on;
for m = 1:size(thres_contrast,1)
    if m == 1
        polyrank = 2;
    else
        polyrank = 3;
    end
    poly_cnt_itr = polyfit(keys_cpd(valid_indx)',10.*1./thres_contrast(m,valid_indx)',polyrank);
    poly_thres_contrast_itr = polyval(poly_cnt_itr,(keys_cpd(1):0.01:keys_cpd(end)));
    h_plt1 = plot(cnt_ax,keys_cpd(valid_indx)',10.*1./thres_contrast(m,valid_indx)','o');
    set(h_plt1,'LineWidth',1,'Marker','o');
    h_plt1.Color(4) = 0.01;
    c = get(h_plt1,'color');
    h_plt2 = plot(cnt_ax,(keys_cpd(1):0.01:keys_cpd(end)), poly_thres_contrast_itr,'--');
    set(h_plt2,'Color',c,'LineWidth',2);
    h_plt1.Color(4) = 0.01;
end
hold off

[~,objects] = legend('CSF all cpd');
objects(2).LineStyle = '-';
objects(2).LineWidth  = 4;
objects(2).Color = 'r';
%% create legend

figure(1);
p1 = plot(ones(100,1),'--o','LineWidth',1,'Marker','o');hold on;
c1 = [0.9294    0.6941    0.1255];
set(p1,'Color',c1,'LineWidth',2);
p2 = plot(2*ones(100,1),'--o','LineWidth',1,'Marker','o');
c2 = [0.4667    0.6745    0.1882];
set(p2,'Color',c2,'LineWidth',2);
p3 = plot(3*ones(100,1),'--o','LineWidth',1,'Marker','o');
c3 = [0    0.4471    0.7412];
set(p3,'Color',c3,'LineWidth',2);
[~,objects] = legend('CSF group 1', 'CSF group 2', 'CSF group 3');
objects(1).FontSize= 18;
objects(2).FontSize= 18;
objects(3).FontSize= 18;

