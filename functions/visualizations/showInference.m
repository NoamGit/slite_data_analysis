function h = showInference(time, f_raw, f, c, s, param)
% initializations
if(nargin < 6)
    lw = 1; % linewidth
else
    if isfield(param,'lineWidth'); lw = param.lineWidth;else lw =1; end
    if isfield(param,'artifact'); arti = param.artifact;else arti = 0; end
    if isfield(param,'baseline'); baseline = param.baseline;else baseline = 0; end
end

    h = figure();
    subplot(311);
    h5 = plot(time,f_raw,'LineWidth',lw);
    axis tight;
    title('F_{raw} (with arti)','FontSize',12);
    set(gca,'FontSize',13,'Box','off');

    subplot(312);
    h3 = plot(time,f,'k','LineWidth',lw);
    title('dF_{processed}','FontSize',12);
    h3.Color(4)=1;
    hold on;
    h4 = plot(time,c,'Color',[1 0.600000023841858 0.7843137383461]...
        ,'LineStyle','--','LineWidth',lw+0.5);
    h4.Color(4)=0.8;
    set(gca,'FontSize',14,'Box','off');
    ylabel('\Delta{F}/F');
    axis('tight')
    hold off;

    subplot(313);
    s(s == 0) = NaN;
    h1 = stem(time,s,'Color',[0.494117647409439 0.494117647409439 0.494117647409439]);
    h1.ShowBaseLine ='off';
    h1.BaseValue = baseline;
    hold on; 
    plot(time,f,'Color',[0.709215695858002 0.909803926944733 0.850980401039124]);
    title('Spike inference with MCMC','FontSize',12);
    h2 = plot(time,c,'Color',[0.584313750267029 0.388235300779343 0.388235300779343],...
        'LineWidth',lw);
    h2.Color(4)=0.8;
   
    set(gca,'FontSize',13,'Box','off');
    xlabel('Time [sec]');
    ylabel('normalized C');
    if(any(arti))
        lw_fill = baseline * ones(size(arti,1),1);
        hg_fill = ones(size(arti,1),1);
        lw_fill(arti <  0.1)= baseline;
        hg_fill(arti <  0.5) = baseline;
        if(size(time,1) > size(time,2))
            h4 = fill( [time' fliplr(time')],  [hg_fill' lw_fill'], 'b');
        else 
            h4 = fill( [time fliplr(time)],  [hg_fill lw_fill], 'b');
        end
        alpha(.09);
        set(h4,'EdgeColor',[1 1 1]);
        set(h4,'EdgeAlpha',0);
        legend('F','est. C','est n_t','stimulus');
        set( legend,'Position',[0.855546467508258 0.432856164761928 0.10220768709019 0.197000685953818] );
        axis([time(1) time(end) (min([f(:) ;c(:)])-eps) max([(-1 * s(:));  f(:) ;c(:)])+eps]);
    else
        legend('F','est. C','est n_t');
        set( legend,'Position',[0.855546467508258 0.432856164761928 0.10220768709019 0.197000685953818] );
        axis('tight')
    end
    hold off

end