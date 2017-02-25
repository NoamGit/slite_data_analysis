function [ output_args ] = sandbox_fn_general( handles,celd  )
        cdata = handles.table_focus(logical(tindx),:); % cell data
        cdata = cdata(2,:);
        
        disp(['SNR is :',num2str(findSNRfromDf(cdata.data.Df, cdata.data.C, cdata.data.S, mean(cdata.data.mcmc_samples.Cb)))]) 

        fig15 = figure(15);clf
        ax15 = axes;
        h1 = plot(celd.time,cdata.data.Df,'Color',[0.559215695858002 0.859803926944733 0.780980401039124]);
        h1.Color(4)=0.99;
        
        hold on;
        h1.Color(4)=0.65;
        h2 = plot(celd.time,cdata.data.C,'Color',[0.584313750267029 0.388235300779343 0.388235300779343],...
        'LineWidth',1);
        h2.Color(4)=0.8;
        cdata.data.S(cdata.data.S <0.06) = nan;
        h0 = stem(celd.time,1.*cdata.data.S ,'Color',[0.3490    0.2000    0.3294],...
            'BaseValue',mean(cdata.data.mcmc_samples.Cb) ,'ShowBaseLine','off');
        
        hold off
        
        set(gcf,'position',[2100 0 871 250]);
        set(gca,'FontSize',13,'Box','on');
        title('Neural Response','FontSize',12);
        xlabel('Time [sec]');
        ylabel('$\Delta$F/F','Interpreter','latex')
%         ylabel('norm $\hat{Ca^+2}$ signal','Interpreter','latex')
        axis_data = get(gca,'YTick');
        axis([celd.time(1) celd.time(end) axis_data(1) axis_data(end)]);
end

