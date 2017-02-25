function [ handles ] = load_fn(handles,hObject)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% resets all axes
cla(handles.axes1)
cla(handles.axes2)
cla(handles.axes5)
cla(handles.axes4)
path = 'C:\Users\noambox\Documents\Sync\';
if(~isdir( path ))
    path = 'C:\Users\Noam\Documents\Sync\';
end

% load cell_loc list
set( handles.detail_TXT,'String','Loading data structure...' );
[fName, pathName] = uigetfile([path, 'Neural data\table_data\*.mat'],'Select Data');
load([pathName, fName]);
data_mode = 'new';
try
    load([pathName,'maps.mat']);
catch 
    disp('working in old dataframe')
    data_mode = 'old';
end
set( handles.detail_TXT,'String','ready!' );
stim_type = textscan(fName,'%s','delimiter','_');
switch stim_type{1}{1}
    case 'ori'
        data = ori_table;
    case 'crp'
        data = crp_table;
    case 'bdn'
        data = bdn_table;
        stim_src_new = load([path,'stimulus files\src\BDN_stimulus_new_5hz.mat']); 
        stim_src_old = load([path,'Sync\stimulus files\src\BDN_stimulus_5hz.mat']);
        handles.bdn_src.new = stim_src_new;
        handles.bdn_src.old = stim_src_old;
    case 'cnt'
        data = cnt_table;
    case 'wgn'
        data = wgn_table;
        stim_src_new = load([path,'stimulus files\src\WGN_stimulus_new.mat']); 
        stim_src_old = load([path,'stimulus files\src\WGN_stimulus.mat']);
        handles.wgn_src.new = stim_src_new.param.out;
        handles.wgn_src.old = stim_src_old.param.out;
    otherwise
        error('stimulus type is not valid')
end

switch data_mode 
    case 'new'
        datel = cellfun(@(c) textscan(c,'%s','delimiter','_'),data.cell_id(:));
        [date, loc] = cellfun(@(c) c{1:2},datel,'UniformOutput' ,false);
        unique_id = cellfun(@(x,y,z) strjoin({x;y;z},'_'),date,loc,data.fname(:),'UniformOutput',false);
        [cell_list, tmap] = unique(unique_id);
        dates = datetime(cellfun(@(c) datestr(c.date), props_map.values(data.fkey(tmap,:)),'UniformOutput', false));
        [~, tmap] = sort(dates,'descend');
        set(handles.cell_location_LST,'Value',1); % must be set if different loads of lists have different lenght
        set(handles.cell_location_LST,'String',cell_list(tmap));% Update handles structure
        handles.prop_map = props_map;
        handles.stim_map = stim_map;    
        
    case 'old'
        unique_id = cellfun(@(x,y) strjoin({x ;y},'_'),data.id(:),data.fname(:),'UniformOutput',false);
        unique_id = cellfun(@(x) textscan(x,'%s','delimiter','-'),unique_id,'UniformOutput',false);
        unique_id = cellfun(@(x) strjoin(x{1}(1:3),'_'),unique_id,'UniformOutput',false);
        set(handles.cell_location_LST,'Value',1); % must be set if different loads of lists have different lenght
        set(handles.cell_location_LST,'String',unique(unique_id));% Update handles structure
end

% plot Df matrix in axes2
df_mat = createDfMat(data.data);
axes(handles.axes4);cla(handles.axes4);
set(gca,'NextPlot','replacechildren');
imagesc(cell2mat(df_mat')');
title('Calcium estimation matrix');xlabel('time[sec]');ylabel('norm Ca^+ signal')
axis tight

% update log
set(handles.fName_TXT,'String',[pathName, fName]);

% update handles
handles.stim_type = stim_type{1}{1};
handles.fName = [pathName, fName];
handles.table = data;
handles.unique_id = unique_id;
handles.data_mode = data_mode;
end

