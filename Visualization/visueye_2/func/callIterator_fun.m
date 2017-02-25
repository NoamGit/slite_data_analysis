function [handles] = callIterator_fun(handles)
% initialze
count = 1;
stopFlag = 0;
h_cpy = handles;
h_cpy.nographics = 1;  
feature_table = table();

% iter short list names
file_iter = Iterator( h_cpy.cell_location_LST.String );
while file_iter.hasNext
    file_str = file_iter.getNext;
    h_cpy = list_short_fun(h_cpy, file_str{:});
    
    % iter long list
    cell_iter = Iterator( h_cpy.cell_LST.String );
    while cell_iter.hasNext()
        cell_str = cell_iter.getNext();
        
        disp(['cell iter - ',num2str(count),' ',file_str{:},' - ',cell_str{:}]);
        
        % get cell and finalize cell data sturcture
        [celd, param, all_signals, all_stim, all_spikes, h_cpy] = lumpingCellData(h_cpy, cell_str);
        try
            celd = buildCellDataStructure(h_cpy, celd, all_signals, all_stim, all_spikes);
        catch e
            disp(e.message);
            % clear all responses with same cell_id from table 
            idx_finished = ismember(h_cpy.table.cell_id, celd.cell_id);
            h_cpy.table(idx_finished,:) = [];
            h_cpy.unique_id(idx_finished,:) = [];
            count = count+1;
            continue;
        end
        
        % process and get features
        try
            feature_cell = eval([h_cpy.funcName,'(celd, h_cpy)']);
            if ~isempty(feature_cell)
                if(istable( feature_cell ))
                    feature_table = [feature_table; feature_cell];
                else
                    temp_table = table();
                    temp_table(celd.cell_id,:) = struct2table(feature_cell);
                    feature_table = [feature_table; temp_table];
                end
            end
        catch e
            disp(e.message);
        end
        
        % clear all responses with same cell_id from table 
        idx_finished = ismember(h_cpy.table.cell_id, celd.cell_id);
        h_cpy.table(idx_finished,:) = [];
        h_cpy.unique_id(idx_finished,:) = [];
        count = count+1;
    end
    
    % update list and iterator    
    file_list = unique(h_cpy.unique_id);
    try
        file_iter = Iterator( file_list );
    catch e
        if(strcmp(e.identifier,'Iterator:notIterable'))
           stopFlag = 1;
        end
    end
    
    if stopFlag == 1; break;
    elseif stopFlag == 0; continue;
    elseif stopFlag == Nan; error('something went wrong!');end;
end

handles.feature_table = feature_table;      
end