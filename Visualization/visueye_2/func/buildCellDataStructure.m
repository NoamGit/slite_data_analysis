function [ cell_data ] = buildCellDataStructure( handles, cell_data, mSig, mStim, mSp )
%rearrages cell_data table

    % plot analysis
    if(strcmp(handles.data_mode,'new'))
        cell_data.props = {handles.prop_map(cell2mat(cell_data.fkey))};
        cell_data.properties = handles.prop_map(cell2mat(cell_data.fkey));
        if(~ismember({'stim'},cell_data.Properties.VariableNames))
            cell_data.stim = {handles.stim_map(cell2mat(cell_data.fkey))};
        end
    end
    cell_data.data.S(isnan(cell_data.data.S)) = 0;
    cell_data.data.S = reshape(cell_data.data.S,1,[]);
    cell_data.multi_data = struct('multi_signal',{mSig}, 'multi_stimulus', {mStim},'multi_spikes',{mSp});

    if(~isempty(handles.samples2remove))
        cell_data = removeData(cell_data, handles.samples2remove);% remove data according to decay rejection
    end
end

