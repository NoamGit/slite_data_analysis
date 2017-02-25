function [stim_out] = calcOffResponseStim(stim)
    [~,peaks_stim]= findpeaks(abs(diff(stim)));
    peaks_stim = [1, peaks_stim+1];
    peaks_pairs = reshape(peaks_stim,2,[]);
    p_diff = diff(peaks_pairs,1);
    p_off = [peaks_pairs(2,:)+1;peaks_pairs(2,:)+1+p_diff];
    stim_out = stim;
    for k = 1:size(p_off,2)
        temp = stim_out(p_off(1,k):p_off(2,k));
        stim_out(p_off(1,k):p_off(2,k)) = stim_out(peaks_pairs(1,k):peaks_pairs(2,k));
        stim_out(peaks_pairs(1,k):peaks_pairs(2,k)) = temp;
    end
end