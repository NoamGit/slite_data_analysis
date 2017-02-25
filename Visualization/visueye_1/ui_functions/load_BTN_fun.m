function [handles] = load_BTN_fun(handles)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% resets all axes
axes(handles.axes2)
cla(gca);
axes(handles.axes3)
cla(gca);
axes(handles.axes4)
cla(gca);
axes(handles.axes5)
cla(gca);

[fName, pathName] = uigetfile('C:\Users\noambox\Documents\Sync\Neural data\*.mat','Select Data');
handles.file = [pathName, fName];
set(handles.file_name,'String',handles.file);
temp = load(handles.file);
if(isfield(temp,'dataframe'))
    handles.df = temp.dataframe; 
else
    handles.df = temp.df; 
end
clear temp;
set(handles.listbox1,'String',handles.df.cell_table.Properties.VariableNames);% Update handles structure

c_mat = varfun(@(x) z1(x.C), handles.df.cell_table);
matplot_examples = (1:handles.df.properties.numCell);
c_mat_array = reshape(table2array(c_mat),handles.df.properties.signalLenght,[]);
s_mat_norm = reshape(table2array(varfun(@(x) z1(x.S), handles.df .cell_table)), handles.df .properties.signalLenght,[]);
s_mat_norm = sparse(s_mat_norm);

axes(handles.axes3);
imagesc(c_mat_array(:,matplot_examples)');
title('Calcium estimation matrix');xlabel('time[sec]');ylabel('norm Ca^+ conc.')
axes(handles.axes4);
spyc(s_mat_norm(:,matplot_examples)','Parula',1);
% daspect([170 2 1])
title('Raster plot');ylabel('cell');xlabel('time[sec]')

% load stimulus
new_flag = handles.df.properties.isnew;
inputsignal = handles.df.cell_table.Cell_29.C;
stim_type = handles.df.properties.stim_type;
pathStim = 'C:\Users\noambox\Documents\Sync\stimulus files\src\';
pathStim_od = 'C:\Users\noambox\Documents\Sync\stimulus files\stim_od';
fs = handles.df.properties.imaging_rate;
numSamples = handles.df.properties.signalLenght;
switch stim_type
    case 'ORI'
        numRep = 3; % for ORI
    case 'BDN'
        numRep = 1;
    case 'WGN'
        numRep = 2;
    case 'CRP'
        numRep = 5;
    case 'CNT' 
        numRep = 3;
    otherwise
        error('!!! - stimulus type is not recognized !!!');
end

% stimulus dictionary .xlsx file
flag_home = 0; % assuming that we are working from lab
try
    stim_dic = readtable('C:\Users\noambox\Documents\Sync\stimulus files\stim_dictionary.xlsx'); 
catch exception % working home
    stim_dic = readtable('C:\Users\Noam\Documents\Sync\stimulus files\stim_dictionary.xlsx'); 
    pathStim = 'C:\Users\Noam\Documents\Sync\stimulus files\src\';
    pathStim_od = 'C:\Users\Noam\Documents\Sync\stimulus files\stim_od';
    flag_home = 1;
end
list_stim = dir(pathStim);
list_stim = {list_stim(~[list_stim.isdir]).name};
ls_tokens = cellfun(@(x) textscan(x,'%s','delimiter',{'_','.mat'}),list_stim);
ls_ind = cellfun(@(x) ismember(stim_type, x),ls_tokens);
ls_ind_new = cellfun(@(x) ismember('new', x),ls_tokens);
[~,indx] = ismember(numSamples,stim_dic.sampels_resp); 
if(indx) % Do we have a match with length in .xlsx table?
    if(~flag_home)
        d_vec = load(stim_dic.path{indx});
    else
        path_lab = stim_dic.path{indx};
        path_tokens = textscan(path_lab,'%s','delimiter',{'\'});
        [~,loc] = ismember('noambox',path_tokens{:});
        path_tokens{1}{3} = 'Noam';
        path_home = strjoin(path_tokens{1},'\');
        d_vec = load(path_home);
    end
    d_vec = d_vec.durationVector;
elseif(~new_flag || ismember(stim_type,{'ORI','CRP'})) % case we don't have a match
    stim_od = load([pathStim_od,'\',stim_type,'_od.mat']); % 1D rep of stimulus
    stim_od = stim_od.stim_od;
    numSamples = numSamples;
    size_duration = numRep * length(stim_od);
    d_vec = diff(linspace(0,(numSamples-1) * 1/fs,size_duration+2))';
    ls_ind = ls_ind & (~ls_ind_new);
else
    stim_od = load([pathStim_od,'\',stim_type,'_od_new.mat']);
    stim_od = stim_od.stim_od;
    size_duration = numRep * length(stim_od);
    d_vec = diff(linspace(0,(numSamples-1) * 1/fs,size_duration+2)');
    ls_ind = ls_ind & ls_ind_new;
end 

if(any(ls_ind) && ~strcmp(stim_type,'CNT')) % subcase 1.2: Do we have the original PTB stimulus file? 
    load([pathStim_od,'\',stim_type,'_od.mat']); 

    % create stimulus vector in the lenght of reponse according to
    % duration vector
%     T = sum(d_vec);
    stim_rs = linspace(0,(numSamples-1)/fs,numSamples)';
    the_map = cumsum(histcounts(stim_rs,cumsum([0 ;d_vec(2:end)])))';
    stim_od_full = repmat(stim_od,1,numRep);
    stim_od_full = stim_od_full(:);
    handles.stimulus = accumarray(the_map,stim_od_full,[numSamples 1],@mean);
    [r, lags] = xcov(handles.stimulus, handles.stimulus(1:ceil(numSamples/numRep))); 
elseif(strcmp(stim_type,'CNT'))
    % TODO: implement for CNT
    load([pathStim_od,'\',stim_type,'_od.mat']); 
    stim_rs = linspace(0,(numSamples-1)/fs,numSamples)';
    the_map = cumsum(histcounts(stim_rs,cumsum([0 ;d_vec(2:end)])))';
    stim_contrs_od_full = repmat(stim_contrs_od,1,numRep);
    stim_contrs_od_full = stim_contrs_od_full(:);
    stim_phase_od_full = repmat(stim_phase_od,1,numRep);
    stim_phase_od_full = stim_phase_od_full(:);
    handles.stimulus = [accumarray(the_map,stim_contrs_od_full,[numSamples 1],@mean)  accumarray(the_map,stim_phase_od_full,[numSamples 1],@mean)];    
    stimUnite = handles.stimulus(:,1) .* handles.stimulus(:,2);
    [r, lags] = xcov(stimUnite, stimUnite(1:ceil(numSamples/numRep))); 
else % we don't have PTB file
    disp('stimulus name is not found in PTB folder')
%     break;
end

% finding the indexes of stimulus repetitions
[~,pksloc] = findpeaks(r(lags>0),'SortStr','descend'); % 
pksloc = sort(pksloc(1:numRep-1));
handles.stimPartition = [[1;pksloc+1] [pksloc;numSamples]]';
% debugging
% stimu = handles.stimulus;
% figure();
% plot([stimu(handles.stimPartition(2,2):handles.stimPartition(2,2)+300-1) stimu(handles.stimPartition(1,2):handles.stimPartition(1,2)+300-1)]);
% plot([stimu(1306:1306+300-1) stimu(2611:2611+300-1)]);


end

