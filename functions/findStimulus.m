function [ df, d_vec ] = findStimulus( df, param)
% findStimulus takes the data frame and extracts the stimulus according to
% the experiment parameter. Same as load_BTN_fun
% crp_flag - a flag for getting low dimention stimulus for CRP
mode = 'default';
crp_flag = 0;
if(nargin>1)
    if isfield(param,'crp_flag'); crp_flag = param.crp_flag;end;
    if isfield(param,'mode'); mode = param.mode;end;
end
if istable(df)
    if(~any(ismember(df.Properties.VariableNames,'properties')))
        df.properties = df;
    end
end

% check were I am 
if( isdir('C:\Users\noambox\Documents\Sync\Neural data'))
    home_lab = 'noambox';
else
    home_lab = 'Noam';
end

% load stimulus
new_flag = df.properties.isnew;
stim_type = df.properties.stim_type;
pathStim = ['C:\Users\',home_lab,'\Documents\Sync\stimulus files\src\'];
pathStim_od = ['C:\Users\',home_lab,'\Documents\Sync\stimulus files\stim_od'];
fs = df.properties.imaging_rate;
numSamples = df.properties.signalLenght;
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
stim_dic = readtable(['C:\Users\',home_lab,'\Documents\Sync\stimulus files\stim_dictionary.xlsx']); 
pathStim = ['C:\Users\',home_lab,'\Documents\Sync\stimulus files\src\'];
pathStim_od = ['C:\Users\',home_lab,'\Documents\Sync\stimulus files\stim_od'];
list_stim = dir(pathStim);
list_stim = {list_stim(~[list_stim.isdir]).name};
ls_tokens = cellfun(@(x) textscan(x,'%s','delimiter',{'_','.mat'}),list_stim);
ls_ind = cellfun(@(x) ismember(stim_type, x),ls_tokens);
ls_ind_new = cellfun(@(x) ismember('new', x),ls_tokens);
[~,indx] = ismember(numSamples,stim_dic.sampels_resp); 

%% loading or creating d_vec and 

if(indx) % Do we have a match with length in .xlsx table?
    if(strcmp(home_lab,'noambox'))
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
    if(crp_flag && strcmp(stim_type,'CRP'))
        load([pathStim_od,'\',stim_type,'_od_low.mat']); % 1D rep of stimulus
        stim_od = reshape(stim_od,[],1);
    else
        load([pathStim_od,'\',stim_type,'_od.mat']); % 1D rep of stimulus
        if( strcmp(stim_type,'CNT') )
            size_duration = numRep * length(stim_contrs_od);
        else
            size_duration = numRep * length(stim_od);
        end
    end
%     stim_od = stim_od.stim_od;
%     numSamples = numSamples;
    d_vec = diff(linspace(0,(numSamples-1) * 1/fs,size_duration+2))';
    ls_ind = ls_ind & (~ls_ind_new);
else
    load([pathStim_od,'\',stim_type,'_od_new.mat']);
    stim_od = reshape(stim_od,[],1);
%     stim_od = stim_od;
    size_duration = numRep * length(stim_od);
    d_vec = diff(linspace(0,(numSamples-1) * 1/fs,size_duration+2)');
    ls_ind = ls_ind & ls_ind_new;
end 

if(strcmp(mode,'duration'))
    df = 0;
    return
end
%% loading stim_od 

if(any(ls_ind) && ~strcmp(stim_type,'CNT')) % subcase 1.2: Do we have the original PTB stimulus file? 
    if(crp_flag && strcmp(stim_type,'CRP'))
        load([pathStim_od,'\',stim_type,'_od_low.mat']); 
    elseif(~strcmp(stim_type,'BDN'))
        load([pathStim_od,'\',stim_type,'_od.mat']); 
    else
        stim_od = zeros(9000,1);
    end
    % create stimulus vector in the lenght of reponse according to
    % duration vector
%     T = sum(d_vec);
    stim_rs = linspace(0,(numSamples-1)/fs,numSamples)';
    the_map = cumsum(histcounts(stim_rs,cumsum([0 ;d_vec(2:end)])))';
    stim_od_full = repmat(stim_od,1,numRep);
    stim_od_full = stim_od_full(:);
    stimulus = accumarray(the_map,stim_od_full,[numSamples 1],@mean);
    [r, lags] = xcov(stimulus, stimulus(1:ceil(numSamples/numRep))); 
elseif(strcmp(stim_type,'CNT'))
    % TODO: implement for CNT
    if(~new_flag) % is old
        load([pathStim_od,'\',stim_type,'_od.mat']); 
    else % is new
        load([pathStim_od,'\',stim_type,'_od_new.mat']); 
    end
    stim_rs = linspace(0,(numSamples-1)/fs,numSamples)';
    the_map = cumsum(histcounts(stim_rs,cumsum([0 ;d_vec(2:end)])))';
    stim_contrs_od_full = repmat(stim_contrs_od,1,numRep);
    stim_contrs_od_full = stim_contrs_od_full(:);
    stim_phase_od_full = repmat(stim_phase_od,1,numRep);
    stim_phase_od_full = stim_phase_od_full(:);
    stimulus = [accumarray(the_map,stim_contrs_od_full,[numSamples 1],@mean)  accumarray(the_map,stim_phase_od_full,[numSamples 1],@mean)];    
    stimUnite = stimulus(:,1) .* stimulus(:,2);
    [r, lags] = xcov(stimUnite, stimUnite(1:ceil(numSamples/numRep))); 
else % we don't have PTB file
    disp('stimulus name is not found in PTB folder')
%     break;
end

% finding the indexes of stimulus repetitions
[~,pksloc] = findpeaks(r(lags>0),'SortStr','descend'); % 
pksloc = sort(pksloc(1:numRep-1));
stimPartition = [[1;pksloc+1] [pksloc;numSamples]]';

df.stimulus = stimulus;
df.stimPartition = stimPartition;
end

