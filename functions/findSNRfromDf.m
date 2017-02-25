%% find SNR of single trace
% Take Raw data and find Peak DF/F0 response over the 
% standard deviation of the signal during one second period before stimulation)
% load data of DF but without Kalmann filtering

function SNR = findSNRfromDf(F, C, S, bl) 
    [peaks, i_p] = sort(S,'descend');
    i_p(i_p > length(F)-3) = [];
    [~,i_p] = findpeaks(C,'MinPeakDistance',10,'SortStr','descend','NPeaks',5);

    %% finding std of 1 sec before stimulus
    J = -stdfilt(flipud(F), ones(21,1));
    [~,maxj] = max(J); % FIXME - not good baseline method
    error('FIXME - baseline estimation samples bad method');
  	bl_std = std(F(maxj-(0:200)));
%     if(size(bl,1) == 1)
%         SNR = mean(F(i_p))/bl_std;  
%     else
    SNR = mean(F(i_p))/bl_std;
%     end
end