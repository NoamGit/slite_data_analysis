function [ cell_data, param, all_signals, all_stim, all_spikes, handles ] = lumpingCellData(handles, cell_str)
% builds the response-matrix. picks all responses of the same cell to the
% specific stimulus in the data, and 
% returns:      cell_data   - a averaged reponse in a cell_data structure (defined by the table) 
%               param       - baseline and handles obj
%               all_signal  - the response matrix aligned
%               all_stim    - all stim repetitions aligned
%               handles     - handles obj

    scale_vec = [0 1];
    
    [~, tindx] = ismember(handles.table_focus.cell_name, cell_str);
    if ~any(ismember( handles.table.Properties.VariableNames,'cell_id'))
        handles.table(:,'cell_id') = cellfun(@(a,b) strcat(a,'_',b), handles.table.id, handles.table.cell_name,'UniformOutput' , false );
        handles.table_focus(:,'cell_id') = cellfun(@(a,b) strcat(a,'_',b), handles.table_focus.id, handles.table_focus.cell_name,'UniformOutput' , false );
    end
    [cell_data, default_props] = allResponsesGivenCellId(handles, tindx);
    
    if(numel(cell_data.data) > 1 && strcmp(handles.data_mode,'new'))
        cell_data.properties = cellfun(@(c) handles.prop_map(c),cell_data.fkey);
        props_table = struct2table(cell_data.properties);
    %     disp(props_table.signalLenght);
        cell_data.stim = cell(numel(cell_data.data),1);
        if(any(props_table.signalLenght-props_table.signalLenght(1))) % case of different signal length
            disp('LIST:: finding duration vectors');
            param.mode = 'duration';
            [~,d_vec] = arrayfun(@(a) findStimulus(cell_data(a,:),param),(1:numel(cell_data.data)),'UniformOutput',false);
            [~,head]= min(cellfun(@sum,d_vec));
            default_props = handles.prop_map(cell2mat(cell_data.fkey(head)));
            cell_data = resampleCellData(head, d_vec,cell_data);
            cell_data.stim(:) = {handles.stim_map(cell2mat(cell_data.fkey(head)))};
        else % same signal length
            cell_data.stim(:) = {cellfun(@(c) handles.stim_map(c),cell_data.fkey(1))};
        end

        dt = 1/default_props.imaging_rate;
        time = (0:dt:( default_props.signalLenght-1 ) * dt);
        if(get(handles.scale_CHBX,'Value'))
            [mean_celd, all_signals, all_stim, all_spikes]  = meanNormedData(cell_data,scale_vec); % TODO: add scale vec
        else
            [mean_celd, all_signals, all_stim,all_spikes]  = meanNormedData(cell_data); % add scale vec
        end;
    %     mean_celd = meanData(celd);
        cell_data((2:end),:) = []; 
        cell_data.data.Df = mean_celd.F; cell_data.data.C = mean_celd.C; cell_data.data.S = mean_celd.S;
        param.baseline = mean_celd.baseline;
        cell_data.time = time;
        param.handles = handles; 
    else% only 1 data recording
        cell_data = cell_data(1,:);
        if(~ismember({'stim'},cell_data.Properties.VariableNames))
            cell_data.stim = {handles.stim_map(cell2mat(cell_data.fkey))};
        else
            cell_data.stim = {struct('stim',cell_data.stim,'partition',cell_data.stim_partition)};
        end
        if(get(handles.scale_CHBX,'Value'))
            [mean_celd, all_signals, all_stim,all_spikes]  = meanNormedData(cell_data,scale_vec); % TODO: add scale vec
        else
            [mean_celd, all_signals, all_stim,all_spikes]  = meanNormedData(cell_data); % add scale vec
        end
        cell_data.data.Df = mean_celd.F; cell_data.data.C = mean_celd.C; cell_data.data.S = mean_celd.S;
        dt = 1/default_props.imaging_rate;
        time = (0:dt:( default_props.signalLenght-1 ) * dt);
        param.baseline = mean(cell_data.data.mcmc_samples.Cb);
        cell_data.time = time;
        param.handles = handles;
        cell_data.data.S(cell_data.data.S < handles.param_SLDR.Value) = 0;
    end
    % auto threshold 
    if( get(handles.athresh_CHBX,'Value') )
        if(~ismember({'stim'},cell_data.Properties.VariableNames))
            cell_data.stim = {handles.stim_map(cell2mat(cell_data.fkey))};
        end
        [ thres ] = findOptimalThreshold(handles,cell_data,'csd');
        set(handles.thresh_TXT,'String',['Spike thresh value: ',num2str(thres)]);
    end
    cell_data.data.S(cell_data.data.S < get(handles.param_SLDR,'Value')) = NaN;
end