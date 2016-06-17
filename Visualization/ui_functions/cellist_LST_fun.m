function [ handles,hObject ] = cellist_LST_fun(handles,hObject)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

contents = cellstr(get(hObject,'String'));
selection = contents{get(hObject,'Value')};
sc = table2array(handles.df.cell_table(:,selection)); % single cell
dt = 1/handles.df.properties.imaging_rate;
time = (0:dt:( handles.df.properties.signalLenght-1 ) * dt);
baseline = mean( sc.mcmc_samples.Cb );
axes(handles.axes2)

sc.S(sc.S == 0) = NaN;
if(handles.st_CHBX.Value)
    h1 = stem(time,(sc.S+baseline),'Color',[0.6667    0.5922    0.7216],...
        'BaseValue',baseline,'ShowBaseLine','off');
    hold on;
else
    h2 = plot(time,sc.C,'Color',[0.584313750267029 0.388235300779343 0.388235300779343],...
        'LineWidth',1);
    set(h2,'PickableParts','all');
    set(h2,'Interruptible','on');
    set(handles.axes2,'PickableParts','all');
    h2.Color(4)=0.8;
    hold on;
end
if(handles.df_CHBX.Value)
    h3 = plot(time,sc.Df,'Color',[0.559215695858002 0.859803926944733 0.780980401039124]);
    h3.Color(4)=0.65;
end
title('Spike inference with MCMC','FontSize',12);
% if(handles.st_CHBX.Value)
    h2 = plot(time,sc.C,'Color',[0.584313750267029 0.388235300779343 0.388235300779343],...
    'LineWidth',1);
    h2.Color(4)=0.8;
% end
set(gca,'FontSize',13,'Box','on','Color',[0.968627452850342 0.968627452850342 0.968627452850342]);
xlabel('Time [sec]');
ylabel('normalized C');
axis tight

% add stimulus inside
if(isfield(handles,'stimulus'))
    stim = handles.stimulus;
    if(strcmp(handles.df.properties.stim_type,'CNT'))
        stim = handles.stimulus(:,1) .* handles.stimulus(:,2);
    end
    stim_scaled = z1(stim) .* 0.3 * (max([(-1 * sc.S(:));  sc.Df(:) ;sc.C(:)])+eps) + 1.1 * (max([(-1 * sc.S(:));  sc.Df(:) ;sc.C(:)])+eps);
    if(handles.stimplt_CHBX.Value) % case 1:  the we plot simulus on the graph
        lw_fill = baseline * ones(size(stim,1),1);
    %     hg_fill = ones(size(stim,1),1);;
        hg_fill = stim/max(stim);
        hg_fill = hg_fill .* (max([(-1 * sc.S(:));  sc.Df(:) ;sc.C(:)])+eps);
        lw_fill(stim == 0)= baseline;
        hg_fill(stim == 0) = baseline;
        if(size(time,1) > size(time,2))
            h4 = fill( [time' fliplr(time')],  [hg_fill' lw_fill'], 'b');
        else 
            h4 = fill( [time fliplr(time)],  [hg_fill; lw_fill], 'b');
        end
        alpha(.2);
        set(h4,'EdgeColor',[1 1 1]);
        set(h4,'EdgeAlpha',0);
        axis([time(1) time(end) (min([sc.Df(:) ;sc.C(:)])-eps) max(stim_scaled+eps)]);
    else % case 2: plot over the graph
        plot(time,stim_scaled,'Color',[0.5686    0.4824    0.4824]);
        axis([time(1) time(end) (min([sc.Df(:) ;sc.C(:)])-eps) max(stim_scaled+eps)]);
    end
    % legend('F','est. C','est n_t','stimulus');
    % set( legend,'Position',[0.855546467508258 0.432856164761928 0.10220768709019 0.197000685953818] );
    % else
    %     legend('F','est. C','est n_t');
    %     set( legend,'Position',[0.855546467508258 0.432856164761928 0.10220768709019 0.197000685953818] );
    %     axis('tight')
    % end
end
    hold off
end

