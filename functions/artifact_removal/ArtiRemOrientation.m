function [ out_signal, out_arti ] = ArtiRemOrientation(signal, arti, param)
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

if nargin < 3 % default
    on_time = 5;
    filt_bs = 11;
    safty = 30;
    filt_thres = 1e3;
    fs = 10;
    maxl_phase = 16;
    corr_thresh = .7;
else
    if isfield(param,'std_filt_binsize'); filt_bs = param.std_filt_binsize;end;
    if isfield(param,'safty'); safty = param.safty; end
    if isfield(param,'on_time'); on_time = param.on_time; end
    if isfield(param,'filt_thres'); filt_thres = param.filt_thres;end;
    if isfield(param,'fs'); fs = param.fs; end
    if isfield(param,'maxl_phase'); maxl_phase = param.maxl_phase; end
    if isfield(param,'corr_thresh'); corr_thresh = param.corr_thresh;else corr_thresh = .7; end
    if isfield(param,'flip_stim'); arti = -arti; end
end

ind = stdfilt(arti,ones(filt_bs,1));
out_signal = signal;
out_arti = arti;
for k = 1:size(ind,2)
    y_iter = signal(:,k);
    a_iter = arti(:,k);
%     lin_ind = find(ind(:,k) > filt_thres);
    [~,locs]=findpeaks(double(ind(:,k) > filt_thres));
    if((~ismember(1,locs)) && ind(1,k) > filt_thres)
        locs = [1; locs]; % first artifact is not recognized by findpeaks
    end
    binsize_ar = on_time*fs + safty;
    loc_mat = repmat(locs', binsize_ar,1);
    range = repmat((-safty/2:(binsize_ar - safty/2-1)),size(loc_mat,2),1)';
    a_ind = loc_mat + range;
    if(any( a_ind(:,1)<1 )) % excluding bounds indexses
        a_ind(:,1) = a_ind(:,1)-a_ind(1,1)+1;
    end
    if( any( a_ind(:,end)> size(y_iter,1)))
        a_ind(:,end) = a_ind(:,end) - ( a_ind(:,end) - size(y_iter,1));
    end
    
    % phase correction
    df_data = num2cell(y_iter(a_ind),1);
    df_arti = num2cell( a_iter(a_ind),1 );
    ind_cell = num2cell(a_ind,1);
    [ df_arti_phac,df_data_phac,ind_phac, best_wind, correction_flag ] = ... % phac ~ phase corrected
        cellfun(@(a, s, ii) correctPhase(a, s, ii, maxl_phase),...
        df_arti,df_data, ind_cell,'UniformOutput',false);
    [~, best_cell] = max(cell2mat( best_wind ));
    [best_sorted,best_ind_sorted] = sort(cell2mat( best_wind ),'descend');
    if(any(best_sorted > corr_thresh))
        window_ind = best_ind_sorted(best_sorted > corr_thresh);
        window_weigths = best_sorted(best_sorted > corr_thresh)/sum(best_sorted(best_sorted > corr_thresh));
        h_v = zeros(2, length(window_ind));
        for jj = 1:length(window_ind)
            b = df_data_phac{window_ind(jj)}; A = df_arti_phac{window_ind(jj)}; 
            A = [A(:) ones(numel(A),1)];

            %find LMS solution for gain factor 'k' 
            h_v(:, jj) = lsqr((A), b(:));
        end
        h = h_v * window_weigths';
    else
        b = df_data_phac{best_ind_sorted(1)}; A = df_arti_phac{best_ind_sorted(1)}; 
        A = [A(:) ones(numel(A),1)];

        %find LMS solution for gain factor 'k' 
        h = lsqr((A), b(:)); 
    end
    
%     % LMS b signals X artifact
%     b = df_data_phac{best_cell}; A = df_arti_phac{best_cell}; 
%     A = [A(:) ones(numel(A),1)];
%     
%     %find LMS solution for gain factor 'k' 
%     h = lsqr((A), b(:));
    art_ridig = A(:,1) * h(1) + h(2); % A_scaled is the scaled 1D signal for further artifact removal
    
    % subtracting phase df_arti_phac,df_data_phac,ind_phac,
    out_arti(:,k) = a_iter * h(1) + h(2);
    out_arti(cell2mat(ind_phac(:)),k) = cell2mat(df_arti_phac(:))* h(1) + h(2);
    out_arti(:,k) = out_arti(:,k) - abs(min(out_arti(:,k)));
    out_signal(:,k) = y_iter - out_arti(:,k);
    
% %     plot
figure(20);
    subplot(311);plot([b(:),art_ridig(:)]); legend('s','fitted a');
    subplot(312);h1 = plot(y_iter );hold on; plot(out_signal(:,k)); h1.Color(4) = 0.6;hold off;axis tight
    subplot(313);plot([y_iter out_signal(:,k)]);
             hold on;plot(-out_arti(:,k)+h(2));hold off;axis tight;

end


    % plot for figures
%     fig1 = figure(1);
%     if(exist(h_v))
%         best_cell = window_ind(end);
%     end
%     plot(signal(a_ind(:,best_cell)));
%     pl_range = (-30:100);
%     b_plot = y_iter(a_ind(1,best_cell)+pl_range);
%     A_plot = a_iter(a_ind(1,best_cell)+pl_range);
%     art_ridig_plot = A_plot(:,1) * h(1) + h(2);
%     t = (a_ind(1,best_cell)+pl_range)/10;
% 
%     plot(gca,t,b_plot); hold on;
%     hf = plot(gca,t,art_ridig_plot,'LineWidth', 2.5);hold off;
%     hf.Color(4) = 0.75;
%     fSize = 22;
%     title('Artifact removal - $s$ vs. $A\cdot\hat{\beta}$','Interpreter','latex','FontSize',fSize)
%     ylabel('Fluo','Interpreter','latex','FontSize',fSize);
%     xlabel('Time [sec]','Interpreter','latex','FontSize',fSize);
%     legend({'$s$','$A\cdot\hat{\beta}$'},'Interpreter','latex','FontSize',fSize);
%     ylim = get(gca,'Ylim');
%     axis([t(1) t(end) ylim(1) ylim(end)]);
% 
%     fig2 = figure(2);clf(fig2);
%     plot(gca,t,b_plot); hold on;
%     hf = plot(gca,t,b_plot - art_ridig_plot +h(2),'LineWidth', 1.5);hold off;
%     hf.Color(4) = 0.75;
%     fSize = 22;
%     title('Artifact removal','Interpreter','latex','FontSize',fSize)
%     ylabel('Fluo','Interpreter','latex','FontSize',fSize);
%     xlabel('Time [sec]','Interpreter','latex','FontSize',fSize);
%     legend({'$s$','$A\cdot\hat{\beta} - s$'},'Interpreter','latex','FontSize',fSize);
%     ylim = get(gca,'Ylim');
%     axis([t(1) t(end) ylim(1) ylim(end)]);

end

function [ a_out, s_out, ind_out, rank, flag ] = correctPhase(arti, s, ind, maxl_phase)
    % find optimal phase with simple xcorr to correct the artifact (source)
    % according to the signal (target)
    % determine - kernel size, binsize
    
    target=norm_nc(s,5);
    source=norm_nc(arti,5);
    
    % compute cross correlation and find lags displacment
    [Cxy,lags]=xcorr(target,source,maxl_phase);
    [~,n] = max(Cxy);
    delta_i = floor(maxl_phase/2) + lags(n) + 1;
        if(delta_i > maxl_phase+1 || delta_i < 1)
            a_out = arti;
            flag = false;
            s_out = s;
            ind_out = ind;
%             a_corrected = arti.*0;
        else
            z = zeros(maxl_phase,1);
            z(delta_i) = 1;
            z = [0; z; 0];
            a_out = conv(arti,z,'valid');
            flag = true;
            s_out = conv(s,z,'valid');
            ind_out = conv(ind,z,'valid') + lags(n);
            
%             figure(9)
%             subplot(311);plot([target source]);
%             subplot(312);plot([target conv(source,z,'same')]);
%             subplot(313);plot(lags,Cxy);
        end
            rank = corrcoef(a_out, s_out);
            rank = rank(2);
    end
   