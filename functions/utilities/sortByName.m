function [ list_out ,map ] = sortByName( list_in, mode)
%sortByName( list_in, mode) 
% list_in - cell array with content to sort

switch mode
    case 'cell_name'
        cell_name_number = cellfun(@(c) textscan(c,'%s','delimiter','_'),list_in);
        cell_name_number = cellfun(@(c) c(2),cell_name_number);
        cell_name_number = cellfun(@(c) str2double(c),cell_name_number);
        [~, map] = sort(cell_name_number); 
        list_out = list_in(map);
        
    case 'cell_id'
        cell_name_number = cellfun(@(c) textscan(c,'%s','delimiter','_'),list_in);
        cell_name_number = cellfun(@(c) c(end),cell_name_number);
        cell_name_number = cellfun(@(c) str2double(c),cell_name_number);
        [~, map] = sort(cell_name_number); 
        list_out = list_in(map);
        
    case 'datetime'
        dates = datetime(cellfun(@(c) datestr(c),list_in),'UniformOutput', false);
        [~, map] = sort(dates,'descend');
        list_out = list_in(map);
        
    otherwise
        disp('soryByCellId:: illeagel mode');
        list_out = list_in;
        map = 0;
end
end

