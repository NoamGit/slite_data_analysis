
% creating stimulus
% note that we use sample [150:(end-1)] sample from the movie. Therefore we
% delete 149 from the beginning and 1 from the end!!!

new_flag = dataframe.properties.isnew;
inputsignal = dataframe.cell_table.Cell_29.C;
stim_type = dataframe.properties.stim_type;
pathStim = 'C:\Users\noambox\Documents\Sync\stimulus files\src\';
pathStim_od = 'C:\Users\noambox\Documents\Sync\stimulus files\stim_od';
fs = dataframe.properties.imaging_rate;
offset = 15;
numRep = 3; % for ORI

% stimulus dictionary .xlsx file
stim_dic = readtable('C:\Users\noambox\Documents\Sync\stimulus files\stim_dictionary.xlsx'); 

list_stim = dir(pathStim);
list_stim = {list_stim(~[list_stim.isdir]).name};
ls_tokens = cellfun(@(x) textscan(x,'%s','delimiter',{'_','.mat'}),list_stim);
ls_ind = cellfun(@(x) ismember(stim_type, x),ls_tokens);
ls_ind_new = cellfun(@(x) ismember('new', x),ls_tokens);
[~,indx] = ismember(length(input),stim_dic.sampels_resp); 
if(indx) % Do we have a match with length in .xlsx table?
    d_vec = load(stim_dic.path{indx});
    d_vec = d_vec.durationVector;
elseif(~new_flag) 
    stim_od = load([pathStim_od,'\',stim_type,'_od.mat']); % 1D rep of stimulus
    stim_od = stim_od.stim_od;
    size_input = length(input);
    size_duration = numRep * length(stim_od);
    d_vec = linspace(0,(size_input-1) * 1/fs,size_duration)';
    ls_ind = ls_ind & (~ls_ind_new);
else
    stim_od = load([pathStim_od,'\',stim_type,'_od_new.mat']);
    stim_od = stim_od.stim_od;
    size_input = length(input);
    size_duration = numRep * length(stim_od);
    d_vec = linspace(0,(size_input-1) * 1/fs,size_duration)';
    ls_ind = ls_ind & ls_ind_new;
end 

% if(new_flag) % subcase 1.1: stimulus is new or old?
%     ls_ind = ls_ind & ls_ind_new;
% else
%     ls_ind = ls_ind & (~ls_ind_new);
% end

if(any(ls_ind)) % subcase 1.2: Do we have the original PTB stimulus file? 
    fname = list_stim(ls_ind);
    load([pathStim,stim_type,'_od.mat']); 

    % create stimulus vector in the lenght of reponse according to
    % duration vector
    T = sum(d_vec);
    stim_rs = linspace(0,(size_input-1)/fs,size_input)';
    the_map = cumsum(histcounts(stim_rs,cumsum([0 ;d_vec(2:end)])))';
    stim_od_full = repmat(stim_od,1,numRep);
    stim_od_full = stim_od_full(:);
    stim_rs = accumarray(the_map,stim_od_full,[size_input 1],@mean);

else % we don't have PTB file
    disp('stimulus name is not found in PTB folder')
%     break;
end


% %          CREATION OF 'STIM'_od!
%         stimulus = load([pathStim,fname{1}]); % make for each stim seperatly
%         stim_od = mean(reshape(stimulus.M,(size(M,1) * size(M,2)),[]),1);
%         stim_od(stim_od > 0 ) = 1;
%         mask = repmat([270, 315, 360, 45, 90 ,135, 180 , 225],450,1);
%         stim_od(1:numel(mask)) = stim_od(1:numel(mask)) .* mask(:)';

% %         CREATION OF CNT_od!
% stim_contrs_od = mean(reshape(stimulus.M,(size(M,1) * size(M,2)),[]),1);
% stim_phase_od = repmat(150*[0.04, 0.06, 0.08, 0.1, 0.2 ,0.4, 0.6 ],600,1);
% stim_phase_od = stim_phase_od(:);
