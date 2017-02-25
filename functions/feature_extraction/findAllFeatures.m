function [ ] = findAllFeatures( files, path )
%finalDataFiltering_fn is the final filtering function before preforming
% meaningfull analysis

global ori_features
global crp_features

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
%% find all features and plot them for ALL cells.
% since we don't have labels we can split to 2 cases - good and bad

default_quantization_size = 16; % for mi and csd calculation
for k = 1:numel(anchor_files)
    
    filepath = anchor_files{k};
    file = load(filepath);
    
    if(isfield(file,'dataframe'))
        df = file.dataframe;
    elseif(isfield(file,'df'))
        df = file.df;
    end
    
    disp(['******* WORKING ON ',filepath,' ****    ']);
    
    % add filename struct
    if(~isfield( df,'sourceFile'))
        df.properties.sourceFile = filepath;
    end
    if(~isfield( df,'ID'))
        [ id,fname ] = findID(filepath, df.properties.date);
        df.ID = id;
    else
        [~,fname] = findID(filepath, df.properties.date);
    end
    
    % find stimulus
    if(~isfield( df,'stimulus'))
        [ df ] = findStimulus(df, 'low_crp');
    end
    
    % calc features for every cell in df
    range_processing = 1:numel(df.cell_table.Properties.VariableNames);
    num_cell = numel(df.cell_table.Properties.VariableNames);
    [datel,filename] = deal( cell(num_cell,1) );
    datel(:) = {df.ID};
    filename(:) = {fname};
    cell_indx = df.cell_table.Properties.VariableNames';
    [ QI, MI, CSD ] = deal( zeros(num_cell,1) );
    df_copy = df;
    for cell_itr = range_processing
        [qi, mi, csd] = calc_features_fn(df_copy, cell_itr, default_quantization_size);
        
        [QI(cell_itr), MI(cell_itr), CSD(cell_itr)] = deal(qi, mi, csd);
%         features_itr(:,cell_itr) = [qi, mi, csd]';
        disp(['finished cell ',num2str(cell_itr)]);
        
        % update cell structure
        var_name = df_copy.cell_table(1,cell_itr).Properties.VariableNames; 
        temp_struct = table2array(df_copy.cell_table(:,var_name));
        temp_struct.QI = qi; 
%         temp_struct.QI_spk = qi_S;
        temp_struct.MI = mi;
%         temp_struct.xcoval = xcov_val;
        temp_struct.CSD = csd;
        
        %% extra for labeling according to critetria
        ori_label(cell_itr) = ((temp_struct.QI > 0.45) & ...
                ((temp_struct.MI > 0.1) &...
                (temp_struct.CSD > 8 ))); 
        
        %%
        df.cell_table(:,var_name) = [];
        df.cell_table(:,var_name) = table(temp_struct,'VariableNames',var_name);
        
    end
%     label = logical(ori_label)';
%     T = table(datel, filename, cell_indx, QI, MI, CSD,label);
    T = table(datel, filename, cell_indx, QI, MI, CSD);

%     T = table(datel, filename, cell_indx);

    switch df.properties.stim_type
        case 'CRP'
            crp_features = [crp_features; T];
        case 'ORI'
            ori_features = [ori_features; T];
    end
    
    cd(path);
    save(filepath,'df');   
end

end

function [id, filename] = findID(filepath, date)
    path_tokens = textscan(filepath,'%s','delimiter','\');
    location = path_tokens{1}{8};
    filename = path_tokens{1}{end};
    id = [datestr(date),'_',location];
end

