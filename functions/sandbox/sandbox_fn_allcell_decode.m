function [  ] = sandbox_fn_allcell_decode( h, cell_data )
% for allcellular rate
% visualize the repetativness of the averaged response of all cells

df = h.table_focus;
temp = df.data;
temp_S = (reshape([temp.S],[],size(df,1)));
temp_S = bsxfun(@rdivide,temp_S,max(temp_S,[],1)); %normalize to max = 1
temp_S(temp_S < h.param_SLDR.Value) = 0;
temp_S(:,logical(prod(isnan(temp_S),1))) = [];
temp_C = (reshape([temp.C],[],size(df,1)));
temp_C = bsxfun(@rdivide,temp_C,max(temp_C,[],1)); %normalize to max = 1
temp_C(:,logical(prod(isnan(temp_S),1))) = [];
S_all = sum(temp_S,2);
C_all = sum(temp_C,2);
cell_data.data.S = reshape(S_all,1,[]);
cell_data.data.C = reshape(C_all,1,[]);

% folding accumulated response
figure(10)
par = cell_data.stim{:}.partition;
c_fold = cell_data.data.C;
s_fold = cell_data.data.S;
c_fold = arrayfun(@(x) c_fold(par(1,x):par(2,x)),(1:size(par,2)),'UniformOutput',false);
s_fold = arrayfun(@(x) s_fold(par(1,x):par(2,x)),(1:size(par,2)),'UniformOutput',false);
subplot(211);plot(c_fold{1}); hold on; plot(c_fold{2});hold off; legend('rep 1','rep 2');title('Ca signal folded reps');xlabel('Time [sec]');
set(gca,'FontSize',14);
subplot(212);plot(s_fold{1});hold on; plot(s_fold{2});hold off;legend('rep 1','rep 2');title('Spike process folded reps');xlabel('Time [sec]');
set(gca,'FontSize',14);
set(gcf,'Position',[1 500 1200 500]);

% fold data for folded prediction
celd_fold = cell_data;
s_fold_mat = zeros(max(cellfun(@length,s_fold)),numel(s_fold));
for k = 1:size(s_fold_mat,2)
    s_fold_mat(linspace(1,par(2,k)+1- par(1,k),par(2,k)+1- par(1,k)), k) = s_fold{k};
end
celd_fold.data.S = reshape(mean(s_fold_mat,2),1,[]);
celd_fold.props{:}.signalLenght = size(s_fold_mat,1);
celd_fold.time = celd_fold.time(1:size(s_fold_mat,1));
celd_fold.data.full_S = cell_data.data.S;

[kernel,time, all_signals,stimulus] = wgnSummary(h,celd_fold);
plotWgn( h,time,kernel,all_signals,stimulus,celd_fold ); 
end

