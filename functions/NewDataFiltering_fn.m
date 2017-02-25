function [ ] = NewDataFiltering_fn( files, path, cell_list)
%NewDataFiltering_fn( files, path ) takes form current location all cells
%with specified cell_id in c2k


global total_informative_cell
% global total_cell_count
% global logger
warning off

%throw away background cells
filenames = {files.name};
tidx = ~cellfun(@isempty,strfind(filenames,'ORI'));
filenames = filenames([find(tidx) find(~tidx)]); % reorder files so ORI is first

% add informative data
disp(['******* Adding data for ', path ,' ****    ']);
info_cell = addCellFromAllFiles( cell_list, filenames, path );

% data statistics
total_informative_cell = total_informative_cell + info_cell;
end

function [cell_count] = addCellFromAllFiles(cell_list, filenames, filepath)
global ori_table
global bdn_table
global crp_table
global wgn_table
global cnt_table
global key
global stim_map   
global props_map 

files = {};
files = cellfun(@(x) [files; filepath,'\',x],filenames);

for k = 1:numel(files)
    filepath = files{k};
    file = load(filepath);
    if(isfield(file,'dataframe'))
        df = file.dataframe;
    elseif(isfield(file,'df'))
        df = file.df;
    end
    
    % remove background cells from cell_list
    if(strcmp(df.properties.stim_type,'ORI'))
        datel = df.ID;
        c2r = ['Cell_',num2str(df.properties.numCell)]; % cell to remove
        remove_idx = ismember(cell_list, [df.ID,'_',c2r]);
        if(any(remove_idx))
            disp('bg cell removed!')
        end
        cell_list(remove_idx) = [];
    end
    
    % find stimulus if not assigned
    if(~isfield( df,'stimulus'))
        try 
            [ df ] = findStimulus(df, 'low_crp');
        catch exception
            disp(exception);
        end 
    end
    
    % build cell_id table and compare
    if(~exist('datel','var'))
        datel = buildDatel(df);
        disp(' ** datel created');
    end
    cell_table = cellfun(@(c) [datel,'_',c],df.cell_table.Properties.VariableNames,'UniformOutput',false)';
    c2k_idx = ismember(cell_table,cell_list); % cell to keep (c2k)
    c2k = cell_table(c2k_idx);
    
    filename = textscan(filepath,'%s','delimiter','\');
    filename = filename{:};
    f_tokens = textscan(filename{end},'%s','delimiter','_');
    f_tokens = f_tokens{:};
    [fname, cell_id,fkey] = deal( cell(numel(c2k),1) );
    fname(:) = {strjoin(f_tokens(1:end-1),'_')};
    cell_id(:) = c2k;
    fkey(:) = {key};
    cell_name = df.cell_table.Properties.VariableNames(c2k_idx)';
    stim_map(key) = struct('stim',df.stimulus,'partition',df.stimPartition);
    props_map(key) = df.properties;
    

    % update informatic cell table
    try
        temp_struct = table2array(df.cell_table(:,cell_name));
    catch exception
        disp(['!!! - in ',datel,' problem with file ',filename(end) ]);
        disp(exception);
    end
    
    data = temp_struct';
    T = table(cell_id,cell_name, fname, fkey, data);
    key = key + 1;
    
    % add to global structure
    switch df.properties.stim_type
        case 'ORI'
            ori_table = [ori_table ; T];
        case 'BDN'
            bdn_table = [bdn_table ; T];         
        case 'WGN'
            wgn_table = [wgn_table ; T];
        case 'CRP'
            crp_table = [crp_table ; T];
        case 'CNT' 
            cnt_table = [cnt_table ; T];
        otherwise
            error('!!! - stimulus type is not recognized !!!');
    end
end

cell_count = numel(c2k);
end

function datel = buildDatel(df)
    loc_indx = strfind(df.properties.sourceFile,'L');
    location = df.properties.sourceFile(loc_indx:loc_indx+1);
    datel = strjoin([{datestr(df.properties.date)},{location}],'_');
end

