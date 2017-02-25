function [handles, hObject] = list_short_fun(handles, hObject)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if ~ischar(hObject)
    contents = cellstr(get(hObject,'String'));
    selection = contents{get(hObject,'Value')};
else
    selection = hObject;
end

% select right table input
switch handles.data_mode
    case 'new'
        [tidx_1] = cellfun(@(c) strcmp(selection,c),handles.unique_id);
        if(~any(tidx_1));return; end;
        cell_list = sortByName(handles.table(tidx_1,:).cell_name, 'cell_name');
        set(handles.cell_LST,'Value',1); % must be set if different loads of lists have different lenght
        set(handles.cell_LST,'String',cell_list);% Update handles structure
        props = handles.prop_map(cell2mat(handles.table(tidx_1,:).fkey(1)));
        % FIXME
%         tidx_file = ismember(handles.table.fname ,handles.table(tidx_1,:).fname);
%         tidx_id = ismember(handles.table.id ,handles.table(tidx_1,:).id);
        tidx_focus = tidx_1;
        tidx = tidx_1;
        stim = handles.stim_map(cell2mat(handles.table(tidx_focus,:).fkey(1))); % get partition repetition vector
        props = handles.prop_map(cell2mat(handles.table(tidx_focus,:).fkey(1))); 
    case 'old'
        [~,tidx_1] = ismember(selection,handles.unique_id);
        if(~any(tidx_1));return; end;
        [~,tidx] = ismember(handles.table.id,handles.table(tidx_1,:).id);
        cell_list = handles.table(logical(tidx),:).cell_name;
        set(handles.cell_LST,'Value',1);
        set(handles.cell_LST,'String',cell_list);% Update handles structure
        props = handles.table(find(logical(tidx),1,'first'),:).props{:};
        tidx_file = ismember(handles.table.fname ,handles.table(tidx_1,:).fname);
        tidx_id = ismember(handles.table.id ,handles.table(tidx_1,:).id);
        tidx_focus = logical(tidx_file) & logical(tidx_id);
        stim = struct('stim', handles.table(tidx_focus,:).stim(1),'partition', handles.table(tidx_focus,:).stim_partition(1)); % get partition repetition vector
        props = handles.table(tidx_focus,:).props(1); 
        props = props{:};
end

if ~ischar(hObject)
    % update log
    if(isfield(props,'sourceFile'))
        txt = sprintf('Stimulus: %s\nDate: %s\nGcamp: %s\nFs: %d\nnumber of cells: %d\nsignal size: %d\nisnew?: %d\nfile: %s',...
                    props.stim_type, datestr(props.date),props.gcampType,props.imaging_rate,props.numCell,props.signalLenght,props.isnew,props.sourceFile);
    else
       txt = sprintf('Stimulus: %s\nDate: %s\nGcamp: %s\nFs: %d\nnumber of cells: %d\nsignal size: %d\nisnew?: %d',...
                props.stim_type, datestr(props.date),props.gcampType,props.imaging_rate,props.numCell,props.signalLenght,props.isnew);
    end
    set( handles.detail_TXT,'String',txt);

    % plot Df matrix of selected file in axes2
    df_mat = createDfMat(handles.table(tidx_focus,:).data);
    axes(handles.axes4);
    cla(handles.axes5,'reset');
    set(gca,'NextPlot','replacechildren');
    imagesc(cell2mat(df_mat')');
    title('Calcium estimation matrix');xlabel('time[sec]');ylabel('norm Ca^+ conc.')
    axis tight
end

% decide if to trimm stimulus repetitions
S_mat = createSpikeMat(handles.table(tidx_focus,:).data);
partition_index = stim.partition;
rep_ignore = filterReps(cell2mat(S_mat'), partition_index, 10 * props.imaging_rate); % winSize  picked to be 10 times the imaging rate
samples2remove = partition_index(:,rep_ignore); % handle data

% update handles 
handles.table_focus = handles.table(logical(tidx),:);
handles.samples2remove = samples2remove;
end

function [rep_ignore] = filterReps(spike_mat, parti, winSize)
% the intuition is that we wnat to check whether  there is a change in the
% firing rate statistics of a big enough window. piecewise stationary
% assumption in winSize duration

% filter unlikely spike
global_thres = 0.01;
spike_mat(spike_mat < global_thres ) = 0;

% divide according to partition
seper_mat = mat2cell(spike_mat, [parti(2,:) + 1 - parti(1,:)], [size(spike_mat,2)]);

% running window on number of spikes
mean_std = @(x) [mean(x), std(x)];
rep_stats = cell2mat(cellfun(@(c) mean_std(runningSpikeCounter(c)) , seper_mat,'UniformOutput',false));

% cut decay according to stats
include_vec = zeros(size(rep_stats,1),1);
include_vec(1) = 1;
for k = 2:size(rep_stats,1)
    reference = mean(rep_stats(logical(include_vec),:),1);
%     if( (rep_stats(k,1) > reference(1) - 0.5*reference(2) ) && (
%     rep_stats(k,1) < reference(1) + 0.5*reference(2) ) ) %falls in
%     statistics
    if( (rep_stats(k,1) + 0.5 * rep_stats(k,1) > reference(1) - 1*reference(2) ))
        include_vec(k) = 1;
    end
end
% disp(include_vec')
rep_ignore = ~include_vec;

    function [stats] = runningSpikeCounter( s )
        nnz_struct = @(block_struct) sum(nonzeros((block_struct.data)));
        stats = blockproc(s,[1,size(s,2)],nnz_struct,'BorderSize',[round(winSize/2),1],...
                'TrimBorder',false,'PadPartialBlocks',true, 'PadMethod','symmetric');
    end

end




