function [ out_signal, out_arti ] = ArtiRemBDN(signal, arti)
% % This fuction takes an Df/f signal and removes the artifact according to
% fluo footages % from the same ROI with fluerecene. The optimal scale and translation 'k'
% is  estimated to obtain the % MLSE solution for the eq: k_hat =
% argmin_k(signal_df - k * artifact)^2 % function returns y = signal_df -
% k_hat * artifact
%     inputs:
%     y - the input signal matrix
%     on_time - Artifact on time display
%     filt_bs - std filter bin size
%     safety_offset - off set in the artifact phase correction
%     filt_thres - filter threshold for artifact detection
%     fs - sampling rate
%     maxl_phase - maxlags for phase xcorr correction

% Segment into parts with artifact values (possibly to use variance
% filtring). Arrange new data and artifact vectors

out_signal = signal;
out_arti = arti;
for k = 1:size(signal,2)
    y_iter = signal(:,k);
    a_iter = arti(:,k);
    
    % take values for the main modal
    perc_values = [.2 .8];
    quantile_vec = quantile(y_iter, perc_values, 1);
    
    % take only < precentile values
    quantile_vec = quantile(y_iter, perc, 1);
    indx_1 = y_iter < quantile_vec;
    
    % take only values where artifact is in 1 std of first modal in the
    % bimodal noise
    d_mu = 12e3; % empirical value of gauss fir to arti_values
    d_std = sqrt(5.85e6);
    indx_2 = (a_iter < d_mu + d_std) & (a_iter > d_mu - d_std);
    
    % exclude outliers in data
    indx_3 = abs(y_iter) < std(y_iter) ;
    indx = indx_2 & indx_1 & indx_3;
    
    A = a_iter(indx);
    A = [A(:) ones(numel(A),1)];
    p = min(y_iter(indx));
    b = y_iter(indx)+ p + 1;
    
    %find LMS solution for gain factor 'k' 
    h = lsqr(A, b(:));
    arti_ridig = a_iter(indx) * h(1) + h(2); % A_scaled is the scaled 1D signal for further artifact removal
    a_iter = a_iter * h(1) + (h(2)-( p + 1 ));
    out_arti(:,k) = a_iter;
    out_signal(:,k) = y_iter - a_iter;
    
% %     plot
    figure(4);plot([b arti_ridig])
    figure(5);subplot(211);plot(y_iter,'r');axis tight;subplot(212);plot( out_signal(:,k));axis tight
    figure(6);subplot(211);plot([y_iter out_signal(:,k)]);axis tight;
             subplot(212);plot([out_arti(:,k)]);axis tight;
end

   