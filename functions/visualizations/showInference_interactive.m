function out_handles = showInference_interactive(time, F, C, S, param)
% showInference_interactive(time, F, C, S, param)
if(nargin < 5)
    lw = 1; % linewidth
else
    if isfield(param,'lineWidth'); lw = param.lineWidth;else lw =1; end
    if isfield(param,'baseline'); baseline = param.baseline;else baseline = 0; end
    if isfield(param,'handles'); handles = param.handles;else handles = 0; end
    if isfield(param,'axis'); ax_vec = param.axis;else ax_vec = 0; end
    if isfield(param,'scale_s'); scale_s_flag = param.scale_s;else scale_s_flag = 1; end
end

% plot order - S, F and last C

S(S == 0) = NaN;
set(gca,'NextPlot','replacechildren');
if(isstruct(handles))
    if(handles.S_CBX.Value) % check S
        if(scale_s_flag)
            temp_s = mapminmax_nc(S,max(C)+0.2 * std(C), baseline);
        else
            temp_s = S+baseline;
        end
        h0 = stem(time, temp_s ,'Color',[0.6667    0.5922    0.7216],...
            'BaseValue',baseline,'ShowBaseLine','off');
        hold on;
    end
    
    if(handles.F_CBX.Value) % check F
        h1 = plot(time,F,'Color',[0.559215695858002 0.859803926944733 0.780980401039124]);
        h1.Color(4)=0.65;
        hold on;
%     axes(handles.axes1);
%     plot(randi(10,100,1));
    end
    
    if(handles.C_CBX.Value) % check C
        h2 = plot(time,C,'Color',[0.584313750267029 0.388235300779343 0.388235300779343],...
        'LineWidth',1);
        h2.Color(4)=0.8;
        hold on;
    end
    
    % reorder
    if(handles.S_CBX.Value)
        uistack(h0,'top');
    end
    if(handles.C_CBX.Value)
        uistack(h2,'top');
    end
    
    if(isfield(param,'stimulus'))
        if(isstruct(param.stimulus))
            stim = param.stimulus.stim;
        else
            stim = param.stimulus;
        end
        stim_type = param.stim_type;
        if(strcmp(stim_type,'CNT'))
            stim = stim(:,1) .* stim(:,2);
        end
        stim_scaled = z1(stim) .* 0.3 * (max([(-1 * S(:));  F(:) ;C(:)])+eps) + 1.1 * (max([(-1 *S(:));  F(:) ;C(:)])+eps);
        if(strcmp(stim_type,'ORI')) % case 1:  the we plot simulus on the graph
            stim = reshape(stim,[],1);
            lw_fill = baseline * ones(size(stim,1),1);
            hg_fill = stim ~= 0;
            hg_fill(hg_fill == 0) = baseline;
            hg_fill = hg_fill .* max([(-1 *S(:));  F(:) ;C(:)]);
            time = reshape(time,[],1);
            lw_fill(stim == 0)= baseline;
            hg_fill(stim == 0) = baseline;
            hg_fill = reshape(hg_fill,[],1);
            lw_fill = reshape(lw_fill,[],1);
            h4 = fill( [time' fliplr(time')],  [hg_fill' lw_fill'], 'b');
            alpha(.2);
            set(h4,'EdgeColor',[1 1 1]);
            set(h4,'EdgeAlpha',0.5);
        else % case 2: plot over the graph
            plot(time,stim_scaled,'Color',[ 0.2039    0.3020    0.4941],'LineWidth',1);
            ax_vec = [time(1) time(end) (min([F(:) ;C(:)])-eps) max(stim_scaled+eps)];
        end
    end
    
    title('Neural Response','FontSize',12);
    set(gca,'FontSize',13,'Box','on','Color',[0.968627452850342 0.968627452850342 0.968627452850342]);
    xlabel('Time [sec]');
    ylabel('norm $\hat{Ca^+2}$ signal','Interpreter','latex');
    if(isfield(param,'axis'))
        axis(ax_vec);
    else
        max_f = max(F);std_f = std(F); min_f = min(F); max_c = max(C);std_c = std(C);min_c = min(C);
        axis([time(1) time(end) min([min_f, min_c])-0.1*min([std_f std_c]) ...
            (max([max_f,max_c]) + 0.1*max([std_f std_c]))]);
    end
    
else
    stem(time,(S+baseline),'Color',[0.6667    0.5922    0.7216],...
        'BaseValue',baseline,'ShowBaseLine','off');
    hold on;
    h2 = plot(time,C,'Color',[0.584313750267029 0.388235300779343 0.388235300779343],...
        'LineWidth',1);
    h2.Color(4)=0.8;
    h3 = plot(time,F,'Color',[0.559215695858002 0.859803926944733 0.780980401039124]);
    h3.Color(4)=0.65;
    
    title('Spike inference with MCMC','FontSize',12);
    h2 = plot(time,C,'Color',[0.584313750267029 0.388235300779343 0.388235300779343],...
    'LineWidth',1);
    h2.Color(4)=0.8;
    set(gca,'FontSize',13,'Box','on','Color',[0.968627452850342 0.968627452850342 0.968627452850342]);
    xlabel('Time [sec]');
    ylabel('norm $\hat{Ca^+2}$ signal','Interpreter','latex')
    axis tight

    % add stimulus inside
    if(isfield(param,'stimulus'))
        stim = param.stimulus;
        stim_type = param.stim_type;
        if(strcmp(stim_type,'CNT'))
            stim = stim(:,1) .* stim(:,2);
        end
        stim_scaled = z1(stim) .* 0.3 * (max([(-1 * S(:));  F(:) ;C(:)])+eps) + 1.1 * (max([(-1 *S(:));  F(:) ;C(:)])+eps);
        if(strcmp(stim_type,'ORI')) % case 1:  the we plot simulus on the graph
            stim_sorted = sort(unique(stim),'ascend');
            stim_sorted = [1  stim_sorted(5:end) stim_sorted(2:4)];
            stim = stim_sorted(stim) - 1;
            lw_fill = baseline * ones(size(stim,2),1);
            hg_fill = stim/max(stim);
            hg_fill = hg_fill .* (max([(-1 * S(:));  F(:) ;C(:)])+eps);
            time = reshape(time,[],1);
            lw_fill(stim == 0)= baseline;
            hg_fill(stim == 0) = baseline;
            hg_fill = reshape(hg_fill,[],1);
            lw_fill = reshape(lw_fill,[],1);
            h4 = fill( [time' fliplr(time')],  [hg_fill' lw_fill'], 'b');
            alpha(.2);
            set(h4,'EdgeColor',[1 1 1]);
            set(h4,'EdgeAlpha',0.5);
            axis([time(1) time(end) (min([F(:) ;C(:)])-eps) max(stim_scaled+eps)]);
        else % case 2: plot over the graph
            plot(time,stim_scaled,'Color',[0.5686    0.4824    0.4824]);
            axis([time(1) time(end) (min([F(:) ;C(:)])-eps) max(stim_scaled+eps)]);
        end
    end
end
hold off;
end