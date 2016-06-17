function [ data_stats ] = dataStatistics(fName, pathName , data_stats )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

pathTokens = textscan( pathName, '%s','delimiter',{' ','_','\','#','  ','Ret'});
nameTokens = textscan(fName,'%s','delimiter','_');

if(~isempty(find(strcmp( 'gcamp6s',pathTokens{1} ),1)))
    data_stats.gcamp6s = data_stats.gcamp6s + 1;
elseif(~isempty(find(strcmp( 'gcamp6f',pathTokens{1} ),1)))
    data_stats.gcamp6f = data_stats.gcamp6f + 1;
end

df = readtable([pathName fName]); % data frame
time = df.Time;
fs = 1/(double(time(2)-time(1)));
off_period = 15*fs;
data_stats.fs = [data_stats.fs ;fs];
stim_list = [{'ORI'};{'BDN'};{'CNT'};{'CRP'};{'WGN'}];
stim_idx = (ismember(nameTokens{1},stim_list));
stim_type = nameTokens{1}{stim_idx};
data_stats.new = data_stats.new + isempty( strfind(pathName,'old') );
num_cell = size(df,2) - 2;
num_samples = size(df,1) - (off_period-1);

% signal length?    % how many cells?
switch stim_type
    case 'ORI'
        data_stats.ori_Len = [data_stats.ori_Len;num_samples];
        data_stats.ori_Cell = data_stats.ori_Cell + num_cell;
        data_stats.total_Cell = data_stats.total_Cell + num_cell;
    
    case 'CNT'
        data_stats.cnt_Len = [data_stats.cnt_Len;num_samples];
        data_stats.cnt_Cell = data_stats.cnt_Cell + num_cell;
        data_stats.total_Cell = data_stats.total_Cell + num_cell;

   case 'CRP'
        data_stats.crp_Len = [data_stats.crp_Len;num_samples];
        data_stats.crp_Cell = data_stats.crp_Cell + num_cell;
        data_stats.total_Cell = data_stats.total_Cell + num_cell;
        
    case 'BDN'
        data_stats.bdn_Len = [data_stats.bdn_Len;num_samples];
        data_stats.bdn_Cell = data_stats.bdn_Cell + num_cell;
        data_stats.total_Cell = data_stats.total_Cell + num_cell;

        
    case 'WGN'
        data_stats.wgn_Len = [data_stats.wgn_Len;num_samples];
        data_stats.wgn_Cell = data_stats.wgn_Cell + num_cell;
        data_stats.total_Cell = data_stats.total_Cell + num_cell;
        
    otherwise
        disp([' **** ',pathName,fName,' has no type?!', char(10) ]);
end

end

