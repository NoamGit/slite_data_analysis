function [ out_sig, out_art ] = ArtiRemFullField(sig, art, param)
% % This fuction takes an Raw signal and removes the artifact according to
% fluo footages from an adjoint 'dark field' region in same image. The optimal scale and translation 'k'
% is  estimated to obtain the 
% MLSE solution for the eq: k_hat =
% argmin_k(signal_df - k * artifact)^2 % function returns y = signal_df -
% k_hat * artifact
%     inputs:
%     sig - the input signal matrix
%     art - the artifact matrix
%     show - plot graphs

if nargin < 3
    show_plot = 0;
    win_size = 200;
else
    if isfield(param,'show_plot'); show_plot = param.show_plot; else show_plot = 0;end
    if isfield(param,'win_size'); win_size = param.win_size; end
end

% constants 
param.Fs_filt = 10;  % Sampling Frequency
param.N       = 5;     % Order
param.Fstop   = 0.5;  % Stopband Frequency
param.Astop   = 80;    % Stopband Attenuation (dB)
    
%% take background artifact and find correlation with art

out_art = repmat(art, [1 size(sig,2)]);
for k = 1:size(sig,2)
    % run moving window for correlation coefficient on the data
    sig_k = detrend_nc( sig(:,k), param );
    [b, A, best_range, best_val] = findCorrRunWin(sig_k, art, win_size, show_plot);
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
        
%         fig2 = figure(2);clf(fig2);
%         t = best_range/ 10;
%         plot(gca,t,b_n * std(art) + mean(art)); hold on;
%         hf = plot(gca,t,a_ls * std(art)+ mean(art),'LineWidth', 2.5);hold off;
%         hf.Color(4) = 0.75;
%         fSize = 22;
%         title('Artifact removal - $s$ vs. $A\cdot\hat{\beta}$','Interpreter','latex','FontSize',fSize)
%         ylabel('Fluo','Interpreter','latex','FontSize',fSize);
%         xlabel('Time [sec]','Interpreter','latex','FontSize',fSize);
%         legend({'$s$','$A\cdot\hat{\beta}$'},'Interpreter','latex','FontSize',fSize);
%         ylim = get(gca,'Ylim');
%         axis([t(1) t(end) ylim(1) ylim(end)]);
%         
%         fig3 = figure(3);clf(fig3);
%         plot(gca,t,b_n* std(art)+ mean(art)); hold on;
%         hf = plot(gca,t,(b_n - A_n(:,1) * h(1)) * std(art)+ mean(art),'LineWidth', 1.5);hold off;
%         hf.Color(4) = 0.75;
%         fSize = 22;
%         title('Artifact removal','Interpreter','latex','FontSize',fSize)
%         ylabel('Fluo','Interpreter','latex','FontSize',fSize);
%         xlabel('Time [sec]','Interpreter','latex','FontSize',fSize);
%         legend({'$s$','$s - A\cdot\hat{\beta}$ '},'Interpreter','latex','FontSize',fSize);
%         ylim = get(gca,'Ylim');
%         axis([t(1) t(end) ylim(1) ylim(end)]);
        end

        out_art(:,k) = (norm_nc(art,6) .* h(1) + h(2)) * std(b);
        out_art(:,k) = out_art(:,k) + abs(min(out_art(:,k)));
end

%% out_art = art .* h(1);
out_sig =  sig - out_art;

% plot
if(show_plot)
    figure();
    for z = 1:3
        subplot(3,1,z)
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