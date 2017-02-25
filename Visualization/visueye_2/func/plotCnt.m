function [out] = plotCnt(handles,mean_response,keys_cpd, varargin )
%UNTITLED4 Summary of this function goes here
    % process response matrix
    proc_response = zeros(size(mean_response));
    keys_cnt = [0 0.1 0.3 0.5 0.7 0.9];
    figure(4);
    thres_contrast = zeros(numel(keys_cpd),1);
    [params, y_hat,f, thres_val] = deal({});
    for k = 1:numel(keys_cpd)
%         curve_to_fit = mean_response((2:end),k);
        imax = 5;
        fixed_params = [0, NaN , NaN , NaN];
        [params{k},~,y_hat{k}, f{k}] = sigm_fit(keys_cnt(2:imax+1)',mean_response((2:imax+1),k),fixed_params,[],0);
        warning off;
        proc_response = f{k}(params{k}(isnan(fixed_params)),keys_cnt);
        subplot(2,4,k);cla;
        if nargin > 5
            errorbar(keys_cnt,mean_response(:,k), varargin{1}(:,k),'LineWidth',2);hold on;
            plot(keys_cnt,proc_response,'--');
        else
            plot(keys_cnt,mean_response(:,k),'o');hold on; plot(linspace(keys_cnt(1),keys_cnt(end),150)...
            ,interp1(keys_cnt,proc_response,linspace(keys_cnt(1),keys_cnt(end),150)),'--');
        end
        [~, ix]= knee_pt(proc_response,keys_cnt',true);
        x50_val = f{k}(params{k}(isnan(fixed_params)),params{k}(3));
        if( params{k}(4) > 1.1 && abs(params{k}(2)) > 5e-3 && (x50_val <= prctile(mean_response(:,k),85)) )
            thres_contrast(k) = params{k}(3);
            thres_val{k} = f{k}(params{k}(isnan(fixed_params)),thres_contrast(k));
        elseif(params{k}(4) > 0.8 && params{k}(4) < 1.1 && abs(params{k}(2)) > 5e-3 && (x50_val <= prctile(mean_response(:,k),85)) )
            [~,ix]= knee_pt(y_hat{k},keys_cnt(2:imax+1)');
            thres_contrast(k) = ( keys_cnt(ix+1)+keys_cnt(ix+2) )/2;
            thres_val{k} = (f{k}(params{k}(isnan(fixed_params)),keys_cnt(ix+1))+f{k}(params{k}(isnan(fixed_params)),keys_cnt(ix+2)))./2;
        elseif(params{k}(2) ~=0 && (x50_val > prctile(mean_response(:,k),85)))
%              thres_contrast(k) = prctile(mean_response(:,k),80);
              [~,ix]= knee_pt(y_hat{k},keys_cnt(2:imax+1)');
              thres_contrast(k) = ( keys_cnt(ix+1)+keys_cnt(ix+2) )/2;
              thres_val{k} = (f{k}(params{k}(isnan(fixed_params)),keys_cnt(ix+1))+f{k}(params{k}(isnan(fixed_params)),keys_cnt(ix+2)))./2;
        elseif(params{k}(2) ==0 || params{k}(4) < 0.5)
            thres_contrast(k) = NaN;
            thres_val{k} = NaN;
        else
            thres_contrast(k) = NaN;
            thres_val{k} = NaN;
        end
        scatter(thres_contrast(k),thres_val{k},'LineWidth',2);hold off
        title(['cpd: ',num2str(keys_cpd(k))])
        xlabel('contrast %');
        ylabel('mean spike values')
    end
    
    figure(5);clf;cnt_ax = axes();
    valid_indx = ~isnan(thres_contrast);
    if sum(valid_indx) >1
    %     spline_thres_contrast = spline(keys_cpd,10.*1./thres_contrast,(min(keys_cpd):0.01:max(keys_cpd)));
        poly_cnt = polyfit(keys_cpd(valid_indx)',10.*1./thres_contrast(valid_indx),2);
        poly_thres_contrast = polyval(poly_cnt,(keys_cpd(1):0.01:keys_cpd(end)));
        h_cnt = loglog(cnt_ax,keys_cpd(valid_indx),10.*1./thres_contrast(valid_indx),'o');
        h_cnt.MarkerSize = 10;
        h_cnt.MarkerFaceColor = [0.75,0.75,0.75];
        grid on
    %     axis([0 0.7 10 50])
        hold on;plot(cnt_ax,(keys_cpd(1):0.01:keys_cpd(end)), poly_thres_contrast); hold off
        title('Contrast sensitivity function')
        xlabel('Spatial freq. [cpd]');
        ylabel('Sensetivity');
        set(gca,'FontSize',20);
    end
    out = struct('params',params,'f',f,'thres_val',thres_val,'thres_contrast',thres_contrast,'y_hat',y_hat);

%   Detailed explanation goes here
if ~isempty(handles)
    axes(handles.axes2);
    set(handles.axes2,'NextPlot','replacechildren');
    imagesc(flipud(mean_response));
    set(handles.axes2,'YTickLabel',fliplr([0 0.1 0.3 0.5 0.7 0.9]))
    set(handles.axes2,'XTickLabel',keys_cpd)
    axis tight;

    axes(handles.axes5);cla;
    set(handles.axes5,'NextPlot','replacechildren');
    copyobj(allchild(cnt_ax),handles.axes5);
    grid on
    title('Contrast sensitivity function')
    xlabel('Spatial freq. [cpd]');
    ylabel('Sensetivity');
    set(gca,'FontSize',12);
    axis([get(cnt_ax,'XLim') get(cnt_ax,'YLim')])
end
end



