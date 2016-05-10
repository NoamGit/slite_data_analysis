%% 150801 Quick illustration script for Shy

source_dir = uigetdir('D:\# Projects (Noam)\# SLITE\# DATA\');
D = dir(fullfile(source_dir, '*.xlsx'));
cd(source_dir);

stimFlag = 0;
cell = {};
for n = 1:length(D)
  if strcmp(D(n).name , 'stimulus.xlsx')
      stim = xlsread(fullfile(source_dir, D(n).name));
      stim = stim(2:end,:);
      stimFlag = 1;
  else
    cell{n} = xlsread(fullfile(source_dir, D(n).name));
  end
end

figure
for n = 1:numel(cell)
    subplot(3,2,n)
    ar = 10; % [hz]if ~isempty(stimulus)
    set(gca,'FontSize',14);
    if stimFlag
        plot(cell{n}(:,1) * 1/ar, (cell{n}(:,2)),'b',cell{n}(:,1) * 1/ar, stim(:,2),'r');
        legend('Ca activity','stimulus','Location','BestOutside');
    else
        plot(cell{n}(:,1) * 1/ar, (cell{n}(:,2)),'b');
    end
        title(['Fluorecence time series - cell ', num2str(n)]); ylabel('\DeltaF/F'); xlabel('time [sec]');
        axis([0 cell{n}(end,1)* 1/ar -inf inf]);
end