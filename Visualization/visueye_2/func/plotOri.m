function [  ] = plotOri(handles, tcd , cell_data )

        % ON tuning curve
        fig3 = figure(3);clf;
        [labels, imap] = sortDegrees(tcd.label);
        [~,where_blank] = ismember('blank',tcd.label);
        [~,where_0] = ismember('0 deg',tcd.label);
        imap(imap==where_blank) = [];
        x = tcd.data(imap,1);
        x = [x; x(1)];labels(end) = {'360 deg'};
        imap =[imap;where_0];
        degs = (0:pi/4:pi*2)';
%         degs = (-pi:pi/4:pi)';
        [~,imax] = max(mean(reshape([x; x; x],3,[]),1));
        tuning_cntr = mod(imax,numel(x));
        if(tuning_cntr == 0);tuning_cntr = 8; end;
        gof.rsquare = 0;
        shift = 0;
        gauss_index = 1; 
        gauss_fit = {'gauss1','gauss2'};
        rand_tuning = 1;
        if(~rand_tuning)
            disp('finding best shift for tuning curve...');
            for k = -5:5 % fixme
                pusher_itr = [3 0 -3 2 -1 -4 1 -2]+k; % TODO: try some offsets
                bm_itr = pusher_itr(tuning_cntr);% bookmark
                x_itr = circshift(x,[bm_itr 0]);
                try
                    [f1,gof1] = fit(degs,x_itr,'gauss1');
                catch
                end
                try
                    [f2,gof2] = fit(degs,x_itr,'gauss2');
                catch
                end
                all_f = {f1,f2};
                all_gof = [gof1, gof2];
                [gof_itr,gauss_index_itr] = max([gof1.rsquare, gof2.rsquare]);
                if(gof.rsquare < gof_itr)
                    gauss_index = gauss_index_itr;
                    gof = all_gof(gauss_index);
                    shift = k;
                    f = all_f{gauss_index};
                end
            end
        else
          gauss_index = 1;
          shift = 0;
        end
        pusher = [3 0 -3 2 -1 -4 1 -2]+shift; % TODO: try some offsets
        bm = pusher(tuning_cntr);% bookmark
        imap_optimal = circshift(imap,[bm 0]);
        x_optimal = circshift(x,[bm 0]);
        [f,gof] = fit(degs,x_optimal,gauss_fit{gauss_index});

        degs_interp = linspace(degs(1),degs(end),200);
        cfit = feval(f,degs_interp);
        std_bar = tcd.data(imap_optimal,2);
        if(min(x_optimal-std_bar) < 0 )
            [~,ind] = min(x_optimal-std_bar);
            htun = errorbar(degs,x_optimal-x_optimal(ind)+std_bar(ind),std_bar,' o');
            hold on;
            plot(degs_interp, cfit-x_optimal(ind)+std_bar(ind));
            hold off;
        else
            htun = errorbar(degs,x_optimal,std_bar,' o');
            hold on;
            plot(degs_interp, cfit);
            hold off;
        end

        set(gca,'XTick',degs,'XTickLabel',[circshift(labels,[bm 0]);{''}],'XTickLabelRotation',90);            
        title('ON response tuning curve','FontSize',14)
        ylabel('Spike counts','FontSize',11);
        set(fig3,'position',[ 1200 500 450 410]);

        % multi signal plot
        fig4 = figure(4);clf;
        data_mat = cell2mat(cell_data.multi_data.multi_signal');
        data_mean = mean(data_mat,2);
        dt = 1/cell_data.properties(:).imaging_rate;
        time = (0:dt:(size(data_mean,1)-1) * dt)';
        data_stim = 2.*any(cell2mat(cell_data.multi_data.multi_stimulus'),2);
        ax_multi = axes;
        hdata = plot(time,data_mat,'Color',[0.75    0.75    0.75]);
        axis_ax = get(ax_multi,'YLim');
        hold on;
        lw_fill = zeros(size(data_stim,1),1);
        hg_fill = data_stim;
        hg_fill = reshape(hg_fill,[],1);
        lw_fill = reshape(lw_fill,[],1);
        hstim = fill( [time' fliplr(time')],  [hg_fill' lw_fill'],'b');
        set(hstim,'EdgeColor',[1 1 1]);
        set(hstim,'FaceColor',[0.8706    0.9216    0.9804]);
        hmean = plot(time,data_mean,'k','LineWidth',1.5);
        axis([time(1) time(end) axis_ax])
        hold off;
        uistack(hdata,'top');
        uistack(hmean,'top');
        xlabel('Time [sec]','Interpreter','latex','FontSize',20);
        ylabel('$\hat{Ca^{+2}}$ signal','Interpreter','latex','FontSize',20);
        title(['0^{\circ}              45^{\circ} ',...
            '            90^{\circ}            135^{\circ}'...
            '           180^{\circ}            225^{\circ}'...
            '           270^{\circ}           315^{\circ}          ']...
            ,'Interpreter','tex')

        set(gca,'FontSize',12);
        set(fig4,'position',[ 1000 200 900 350]);
end

function [y, map] = sortDegrees(x)
    x_num = zeros(numel(x),1);
    for c = 1:numel(x)
        x_scan = textscan(x{c},'%s','delimiter',' ');
        x_num(c) = str2double(cell2mat(x_scan{1}(1)));
    end
    [val,map] = sort(x_num);
    y = cellfun(@(c) strjoin([{num2str(c)},{'deg'}],' '), num2cell(val(1:end-1)),'UniformOutput',false);
    y = [y;{'blank'}];
end