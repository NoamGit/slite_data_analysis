function [  ] = sandbox_fn_foldCnt( h, cell_data )

% folding accumulated response
figure(10)
par = cell_data.stim{:}.partition;
stim_cpd = norm_nc(cell_data.stim{:}.stim(:,2),1);
stim_contrast = norm_nc(cell_data.stim{:}.stim(:,1),1);
c_fold = cell_data.data.C;
s_fold = cell_data.data.S;
c_fold = arrayfun(@(x) c_fold(par(1,x):par(2,x)),(1:size(par,2)),'UniformOutput',false);
stim_cpd_fold = arrayfun(@(x) stim_cpd(par(1,x):par(2,x)),(1:size(par,2)),'UniformOutput',false);
stim_cnt_fold = arrayfun(@(x) stim_contrast(par(1,x):par(2,x)),(1:size(par,2)),'UniformOutput',false);

hg_fill = stim_cpd_fold{1};
lw_fill = zeros(size(hg_fill));
s_fold = arrayfun(@(x) s_fold(par(1,x):par(2,x)),(1:size(par,2)),'UniformOutput',false);
mSig = cell2mat(cell_data.multi_data.multi_signal');
subplot(211);cla;plot(mSig); hold on; 
plot(mean(mSig,2),'LineWidth',2,'Color','k');
h_stim = fill( [(1:length(lw_fill)) fliplr(1:length(lw_fill)')],  [max(mSig(:)) .* hg_fill' lw_fill'], 'b');
alpha(.2);
set(h_stim,'EdgeColor',[1 1 1]);
set(h_stim,'EdgeAlpha',0.5);

h_stim2 = fill( [(1:length(lw_fill)) fliplr(1:length(lw_fill)')],  [max(mSig(:)) .* stim_cnt_fold{1}' lw_fill'], 'r');
alpha(.15);
set(h_stim2,'EdgeColor',[1 1 1]);
set(h_stim2,'EdgeAlpha',0.15);

hold off; legend('rep 1','rep 2','rep 3');title('Ca signal folded reps');xlabel('Time [sec]');
set(gca,'FontSize',14);

subplot(212);cla;plot(s_fold{1});hold on;
plot(s_fold{2});
plot(s_fold{3});
h_stim = fill( [(1:length(lw_fill)) fliplr(1:length(lw_fill)')],  [max(cellfun(@max, s_fold)).*hg_fill' lw_fill'], 'b');
alpha(.2);
set(h_stim,'EdgeColor',[1 1 1]);
set(h_stim,'EdgeAlpha',0.5);

h_stim2 = fill( [(1:length(lw_fill)) fliplr(1:length(lw_fill)')],  [max(cellfun(@max,s_fold)) .* stim_cnt_fold{1}' lw_fill'], 'r');
alpha(.15);
set(h_stim2,'EdgeColor',[1 1 1]);
set(h_stim2,'EdgeAlpha',0.15);

hold off;legend('rep 1','rep 2','rep 3');title('Spike process folded reps');xlabel('Time [sec]');
set(gca,'FontSize',14);
set(gcf,'Position',[1 500 1200 500]);

end
