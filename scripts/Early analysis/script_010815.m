%% 150801 Quick illustration script for Shy

source_dir = uigetdir('D:\# Projects (Noam)\# SLITE\# DATA\');
D = dir(fullfile(source_dir, '*.txt'));
cd(source_dir);

stimFlag = 0;
cell = {};
for n = 1:length(D)
  if strcmp(D(n).name , 'stimulus.txt')
      txtdata = importdata(fullfile(source_dir, D(n).name));
      stim = txtdata.data(2:end,2);
      stimFlag = 1;
  else
    txtdata = importdata(fullfile(source_dir, D(n).name));
    samples = txtdata.data(:,1);
    cell{n} = txtdata.data(:,2);
  end
end

figure
for n = 1:numel(cell)
    subplot(3,2,n)
    ar = 10; % [hz]if ~isempty(stimulus)
    set(gca,'FontSize',14);
    if stimFlag
       if mean(cell{n})>10
            plot(samples* 1/ar, (cell{n}),'b',samples * 1/ar, stim * max(cell{n}),'r');
       else    
            plot(samples* 1/ar, (cell{n}),'b',samples * 1/ar, stim,'r');
       end
        legend('Ca activity','stimulus','Location','BestOutside');
    else
        plot(samples * 1/ar, (cell{n}),'b');
    end
        title(['Fluorecence time series - cell ', num2str(n)]); ylabel('\DeltaF/F'); xlabel('time [sec]');
        axis([0 samples(end)* 1/ar -inf inf]);
end