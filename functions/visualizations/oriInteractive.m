function txt = oriInteractive(empt,event_obj)
% Customizes text of data tips
global ori_features_filt;
global lbl;
global ori_table;
global ori_features;

pos = get(event_obj,'Position');
[~, indx] = ismember( [pos(1),pos(2),pos(3) ],[ ori_features_filt.QI,ori_features_filt.MI,ori_features_filt.CSD ],'rows');
location_txt = ['QI: ',num2str(pos(1)),'  MI: ',num2str(pos(2)),'  CSD: ',num2str(pos(3))];
txt = {['File: ',cell2mat( lbl(indx,1) )],['Location: ',cell2mat( lbl(indx,2))],location_txt};
cell_indx = ori_features_filt.cell_indx(indx);
disp(['cell index is - ',cell_indx]);
[ loc1 ] = ismember( ori_table.id, lbl(indx,1) );
[ loc2 ] = ismember(  ori_table(loc1, :).fname,lbl(indx,2) );
temp_table = ori_table(loc1, :);
data_indx = strcmp(temp_table(loc2, :).cell_name,cell_indx);
all_indx = strcmp(temp_table.cell_name,cell_indx);

all_cells = temp_table(all_indx,:);
temp_table = temp_table(loc2,:);
df = temp_table(data_indx,:);
df_data = temp_table(data_indx,:).data;
l = temp_table(data_indx,:).props{1}.signalLenght;
dt = 1/temp_table(data_indx,:).props{1}.imaging_rate;
time = linspace(0,(l-1) * dt,l);

num_val_stim = 9;
[partition, ~] = lloyds(df.stim{:},num_val_stim);
[~,stim] = quantiz(df.stim{:},partition,(1:num_val_stim));
[all_data_S, all_data_C, all_data_F] = deal(zeros(1,l));
all_baseline = 0;
num_reps = numel(all_cells(:,1));
for k = 1:num_reps
%     disp([num2str(numel(all_cells(:,1))),' reps concluded...']);
    if(length( all_cells(k,:).data.S) ~= length(all_data_S))
        disp('CAUTION! - reps have different size');
        num_reps = k-1;
        break
    else
        all_data_C = all_data_C + all_cells(k,:).data.C;
        all_data_S = all_data_S + all_cells(k,:).data.S;
        all_data_F = all_data_F + all_cells(k,:).data.Df';
        all_baseline = all_baseline + mean(all_cells(k,:).data.mcmc_samples.Cb);
    end
end
mean_s = all_data_S./num_reps;
mean_c = all_data_C./num_reps;
mean_f = all_data_F./num_reps;
all_baseline = all_baseline./num_reps;
param.baseline = all_baseline;
polar_vec = accumarray( stim',mean_s,[],@mean );
% S_filt = df_data.S';
% S_filt( S_filt < 0.05 ) = 0;
% polar_vec = accumarray( stim',df_data.C,[],@mean );
% disp(polar_vec');
polar_vec = [polar_vec(7:end)' polar_vec(2:6)'];
polar_vec = polar_vec./max(polar_vec);
% disp(polar_vec);

% polar plot
theta = linspace(0, 2*pi, 9);
rho = [polar_vec';polar_vec(1)];
figure(11)
p = polar(theta,rho');
all_lines = findall(gca,'type','line'); 
polarticks(8,all_lines(1), gca);

% muliple comparison 1-way Anova test
degree_map = {'blank','135 deg','180 deg','225 deg','270 deg','315 deg','0 deg','45 deg','90 deg'};
deg_vec = degree_map(stim);
[p,t,stats] = anova1(mean_s,deg_vec,'off');
figure(12)
[c,m,hh,nms] = multcompare(stats);

param.stimulus = stim;
param.stim_type = df.props{:}.stim_type;
showInference_interactive(time, mean_f, mean_c, mean_s, param)

