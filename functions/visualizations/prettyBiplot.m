function prettyBiplot( handle )

    dots_list= findobj(handle,'Tag','obsmarker');
    for k = 1:numel(dots_list)
        set(dots_list(k),'Marker','o','MarkerSize',4,'MarkerFaceColor',[0.8510    0.3255    0.0980]);
    end
    line_list= findobj(handle,'Tag','varline');
    for k = 1:numel(line_list)
        set(line_list(k),'LineWidth',1.5,'Color',[ 0    0.4471    0.7412],'LineStyle','-.');
    end
    line_marker_list = findobj(handle,'Tag','varmarker');
    for k = 1:numel(line_marker_list)
        set(line_marker_list(k),'Marker','V','MarkerFaceColor',[ 0    0.4471    0.7412],'Color',[ 0    0.4471    0.7412]);
    end
    text_list = findobj(handle,'Type','Text');
    for k = 1:numel(text_list)
        set(text_list(k),'FontSize',13);
    end
end