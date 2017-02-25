set(gca,'FontSize',25);
title('OFF response tuning curve','FontSize',25);ylabel('');
xtic = get(gca,'XTickLabel');
xtic2 = cellfun(@(c) c{1},cellfun(@(c) textscan(c,'%s','delimiter',' ')...
    ,xtic(1:end-1)),'UniformOutput',false);
set(gca,'XTickLabel',xtic2)
%%
savemultfigs

%% fix labels
xtic = get(gca,'XTickLabel');
xtic = [num2str( str2double(xtic(1)) - 45) ; xtic];
set(gca,'XTickLabel',xtic2)



