function [ out_sig, out_art ] = ArtiRemContrast(sig, art, param)
% % This fuction takes an Raw signal and removes the artifact according to
% fluo footages % from the same ROI with fluerecene. The optimal scale and translation 'k'
% is  estimated to obtain the % MLSE solution for the eq: k_hat =
% argmin_k(signal_df - k * artifact)^2 % function returns y = signal_df -
% k_hat * artifact
%     inputs:
%     sig - the input signal matrix
%     art - the artifact matrix
%     safty - for phase correction. maximum lag of signal
%     show - plot graphs

if nargin < 3
    show_plot = 0;
    safty = 5;
    win_size = 200;
else
    if isfield(param,'show_plot'); show_plot = param.show_plot; else show_plot = 0;end
    if isfield(param,'safty'); safty = param.safty; end
    if isfield(param,'win_size'); win_size = param.win_size; end
end

% constants 
lev = 3; % wavlets denoising level

% take background artifact and find correlation with art
art = bsxfun(@plus, art, abs(min(art)));
a = art(:,end);
s = sig(:,end);
[xc, l] = xcorr(s,a,safty,'unbiased');
ind = find(xc == max(xc));
lag = l(ind);
z = zeros(safty*2+1,1);
z(ind) = 1;
z = [0; z; 0];

if(lag ~= 0)
    disp(['artifact and signal are not fully sync, lag =', num2str(lag)]);
    for k = 1:size(art,2)   % alinged art
        art(:,k) = conv(art(:,k),z,'same');
    end
end

out_art = art;
for k = 1:size(sig,2)
    % run moving window for correlation coefficient on the data
    sig_k = wden(sig(:,k),'sqtwolog','s','sln',lev,'sym8');
    [b, A, range] = findCorrRunWin(sig(:,k), art(:,k), win_size, show_plot);
    A_n = norm_nc(A,6);
    b_n = norm_nc(b,6);

    % sum every binn and take the 5 smallest sums that would represent the
    A_n = [A_n(:) ones(numel(A_n),1)];

    %find LMS solution for gain factor 'k' 
    h = lsqr(A_n, b_n(:));
    a_ls = A_n(:,1) * h(1) + h(2); % A_scaled is the scaled 1D signal for further artifact removal

    % some plotting
    if(show_plot)
        figure()
        subplot(211)
        plot(b_n * std(b));hold on; h1 =  plot(a_ls * std(b));
        h1.Color(4) = 0.4;
        hold off
        xlabel('time [sec]'); ylabel('norm F_t'); legend('s','fitted artifact');
        subplot(212)
        h2 = plot( b);hold on; plot(b- a_ls * std(b));
        h2.Color(4) = 0.8;
        hold off;
        xlabel('time [sec]'); ylabel('F_t'); legend('s','s - a');
    end
    
    out_art(:,k) = (norm_nc(art(:,k),6) .* h(1) + h(2)) * std(b);
    out_art(:,k) = out_art(:,k) + abs(min(out_art(:,k)));
end

% out_art = art .* h(1);
out_sig =  sig - out_art;

% plot
if(show_plot)
    figure();
    for z = 1:4
        subplot(4,1,z)
        plot(out_sig(:,z));hold on; 
        h1 =  plot(sig(:,z),'r'); h1.Color(4) = 0.3; 
%         h2 = plot(out_art(:,z),'k');h2.Color(4) = 0.6;
        legend('sig - art', 'original sig');
        hold off
        ylabel('raw F data - artifact');
        axis tight
    end
end

end

function [y1, y2, best_range_rw, best_val_rw] = findCorrRunWin(x1, x2, win_size,show_plot)
    best_range_rw = []; % running window index
    best_val_rw = 0; % best value of correltation coeff in window
    for m = 1:(size(x1,1) - win_size)
        range = (m:m+win_size-1);
        val_rw = corrcoef(x1(range), x2(range));
        if(val_rw(2) > best_val_rw)
           best_range_rw = range; 
           best_val_rw = val_rw(2);
        end
    end
    y1 = x1(range);
    y2 = x2(range);
    
    if(show_plot)
        figure()
        subplot(2,1,1)
        plot(z1(y1));hold on; 
        h1 =  plot(z1(y2)); h1.Color(4) = 0.4;
        hold off;
    end
end