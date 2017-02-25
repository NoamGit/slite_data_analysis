% finalDataFiltering_script
% also is used to generate ori_feattures and crp_features
global ori_features
global crp_features
global loc_ori
global loc_crp
global ori_table
global bdn_table
global crp_table
global wgn_table
global cnt_table
global total_informative_cell
global total_cell_count
global logger
global key
global stim_map   
global props_map   

if( isdir('C:\Users\noambox\Documents\Sync\Neural data'))
    home_lab = 'noambox';
else
    home_lab = 'Noam';
end
load(['C:\Users\',home_lab,'\Documents\Sync\Neural data\stat_features\all_feat_spectral_last.mat'])

stim_map = containers.Map('KeyType','int32','ValueType','any');
props_map = containers.Map('KeyType','int32','ValueType','any');
key = 1;
logger = [];
[ total_informative_cell, total_cell_count ] = deal(0);
[ ori_table,bdn_table,crp_table, wgn_table, cnt_table ]  = deal( table());

loc_ori = 1;
loc_crp = 1;
%% hard thresholding on features ORI yield 224 cells CRP yield 81 cells

ori_label = ((table2array( ori_features(:,4) ) > 0.45) & ...
        ((ori_features(:,5).MI > 0.1) &...
        (ori_features(:,6).CSD > 8 ))); 
ori_features.label = ori_label;
% nnz(ori_label)

crp_label = (( crp_features(:,4).QI  > 0.35) & ...
        ((crp_features(:,5).MI > 0.08) &...
        (crp_features(:,6).CSD > 4.5 )));   
crp_features.label = crp_label;

% nnz(crp_label)
%% turning datel to cell_id

ori_features.datel = cellfun(@(a,b) strjoin([{a};{b}],'_'),ori_features.datel, ori_features.cell_indx,'UniformOutput',false);
crp_features.datel = cellfun(@(a,b) strjoin([{a};{b}],'_'),crp_features.datel, crp_features.cell_indx,'UniformOutput',false);
cells_2_keep = unique([ori_features.datel(logical(ori_features.label)); crp_features.datel(logical(crp_features.label))]); 
f2itr_list = cell(numel(cells_2_keep),1); % files to iterate on 
for k = 1:numel(cells_2_keep)
    cid = cells_2_keep(k);
    cid_token = textscan(cid{:},'%s','delimiter','_');
    f2itr_list{k} = strjoin(cid_token{:}(1:2),'_');
end
%%  Iterator
  
list_parent_path = dir(['C:\Users\',home_lab,'\Documents\Sync\Neural data']);
for parent_path = fliplr({list_parent_path.name})
    if(strfind(parent_path{1}, 'gcamp'))
        list_child_path = dir(['C:\Users\',home_lab ,'\Documents\Sync\Neural data\',parent_path{1}]);
        for child_path = {list_child_path.name}
            if(strfind(child_path{1}, 'L'))
                pathName = ['C:\Users\',home_lab ,'\Documents\Sync\Neural data\',parent_path{1},'\',child_path{1},'\Parsing'];
                fName_table = dir([pathName,'\*.mat']);
                fName = fliplr({fName_table.name});
                
                % for infor/noninfor segmentation
                findAllFeatures(fName_table, pathName); 
%                 finalDataFiltering_fn(fName_table, pathName);    
                
%                 fname_token = textscan(fName{1},'%s','delimiter',{'_','.mat'});
%                 datel = strjoin([fname_token{:}(end);child_path],'_');
%                 if(any(ismember(f2itr_list,datel)))
%                     NewDataFiltering_fn(fName_table, pathName, cells_2_keep);         
%                 end
            end
        end
    end
end
disp([char(10,'all done!')])