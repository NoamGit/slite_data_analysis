function [ ] = finalDataFiltering_fn( files, path )
%finalDataFiltering_fn is the final filtering function before preforming
% meaningfull analysis

global ori_features
global crp_features
global loc_ori
global loc_crp
global total_informative_cell
global total_cell_count
global logger
warning off

% finds target files for filtering (CRP and ORI)
filenames = {files.name};
anchor_files = {};
ori_indx = find(~cellfun(@isempty,(strfind(filenames , 'ORI'))));
if(~isempty(ori_indx))
    ori_file = filenames(ori_indx);
    anchor_files = cellfun(@(x) [anchor_files; path,'\',x],ori_file');
end

anchor_files_crp = {};
crp_indx = find(~cellfun(@isempty,(strfind(filenames , 'CRP'))));
if(~isempty(crp_indx))
    crp_file = filenames(crp_indx);
    anchor_files_crp = cellfun(@(x) [ anchor_files_crp; path,'\',x ],crp_file');
end
anchor_files = [anchor_files;anchor_files_crp];
%% filters both returning 'garbage cells' - Done with Spectral clustering in IPython
% thres right value is for 'ORI' and left value for 'CRP'

labels = cell(numel(anchor_files),1);
for k = 1:numel(anchor_files)
    
    filepath = anchor_files{k};
    file = load(filepath);
    if(isfield(file,'dataframe'))
        df = file.dataframe;
    elseif(isfield(file,'df'))
        df = file.df;
    end

    var_names = df.cell_table.Properties.VariableNames;
    switch df.properties.stim_type
        case 'ORI'
            labels{k} = var_names(logical(table2array( ori_features...
                (loc_ori:loc_ori+df.properties.numCell-1,7))));
            loc_ori = loc_ori + df.properties.numCell;
            cell_counter = df.properties.numCell-1;
            
            % throw highest value cell from labels
            if(~isempty(labels{k})) % if no single good cell in file
                try
                    indx_bg = hasBackground(labels{k},df.properties.numCell);
                catch exception
                   disp('!!! PROBLEM !!!');
                   disp(exception); 
                end
                if(indx_bg)
                   labels{k}(indx_bg) = [];
                   logger = [logger; {filepath}];
                end
            end
        case 'CRP'
            labels{k} =  var_names(logical(table2array(crp_features...
                (loc_crp:loc_crp+df.properties.numCell-1,7))));
            loc_crp = loc_crp + df.properties.numCell;
            cell_counter = df.properties.numCell;
    end       
end

labels = unique([labels{:}]);

% recall - we work in every Parsing library
% add informative data
datel_id = df.ID;
disp(['******* Adding data for ',datel_id,' ****    ']);
addCellFromAllFiles( labels, filenames, path, datel_id );

% data statistics
total_cell_count = total_cell_count + cell_counter ;
total_informative_cell = total_informative_cell + numel(labels);

end

function [] = addCellFromAllFiles(labels, filenames, filepath, datel_id)
global ori_table
global bdn_table
global crp_table
global wgn_table
global cnt_table
global key
global stim_map   
global props_map 

files = {};
files = cellfun(@(x) [files; filepath,'\',x],fliplr(filenames));

for k = 1:numel(files)
    filepath = files{k};
    file = load(filepath);
    if(isfield(file,'dataframe'))
        df = file.dataframe;
    elseif(isfield(file,'df'))
        df = file.df;
    end
    
    % find stimulus if not assigned
    if(~isfield( df,'stimulus'))
        try 
            [ df ] = findStimulus(df, 'low_crp');
        catch exception
            disp(exception);
        end 
    end
    

    % TODO - check if cell's location are simmilar for different files. if
    % not change indexing
    filename = textscan(filepath,'%s','delimiter','\');
    filename = filename{:};
    f_tokens = textscan(filename{end},'%s','delimiter','_');
    f_tokens = f_tokens{:};
    [fname, cell_id, props,stim, stim_partition,fkey] = deal( cell(numel(labels),1) );
    fname(:) = {strjoin(f_tokens(1:end-1),'_')};
    cell_id(:) = cellfun(@(c) strjoin([{datel_id};{c}],'_'), labels,'UniformOutput',false)';
    fkey(:) = {key};
    cell_name = labels';
    stim_map(key) = struct('stim',df.stimulus,'partition',df.stimPartition);
    props_map(key) = df.properties;
    
%     props(:) = {df.properties};
%     stim(:) = {df.stimulus};
%     stim_partition(:) = {df.stimPartition};
        % update informatic cell table
    try
        temp_struct = table2array(df.cell_table(:,labels));
    catch exception
        disp(['!!! - in ',datel_id,' problem with file ',filename(end) ]);
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


end

function [out] = hasBackground(cellA, b)
    cellA = cellfun(@(x) (textscan(x,'%s','delimiter',{'_'})),cellA');
    cellA = [cellA{:}];
    numA = str2double(cellA(2,:));
    [~,out] = ismember(b, numA);
end
