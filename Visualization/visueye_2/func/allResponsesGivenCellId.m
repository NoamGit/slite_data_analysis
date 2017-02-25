function [ cell_data, default_props, all_cell_ids ] = allResponsesGivenCellId(handles, cell_index)

    % extract data
    cell_data = handles.table_focus(logical(cell_index),:); % cell data
    if(strcmp(handles.data_mode,'new'))
        default_props = handles.prop_map(cell2mat(cell_data.fkey));
        tidx = ismember(handles.table.cell_id,cell_data.cell_id);
        cell_data = handles.table(tidx,:);
    else
        default_props = cell_data.props{1};
    end
    
    all_cell_ids = cell_data(:,'cell_id');
    % clean empty traces from data
    for k = fliplr(1:numel(cell_data.data))
        if(~any(cell_data.data(k).S) && numel(cell_data.data) > 1)
            cell_data(k,:) = [];
        end
    end
end