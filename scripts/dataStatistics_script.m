% parse data script 
% moves across all files from a .xlmx type in the directory folder and subfolder and execute the parse
% data function
global logvar;
logvar = 'List of files with no artifact:';

data_stats = struct('gcamp6f',0,'gcamp6s',0,'fs',[],'new',0,'ori_Len',[],'ori_Cell',0,...
                    'cnt_Len',[],'cnt_Cell',0,'crp_Len',[],'crp_Cell',0,'bdn_Len',[],...
                    'bdn_Cell',0,'wgn_Len',[], 'wgn_Cell',0,'total_Cell',0);

warning off;
list_parent_path = dir('C:\Users\noambox\Documents\Sync\Neural data');
for parent_path = fliplr({list_parent_path.name})
    if(strfind(parent_path{1}, 'gcamp'))
        list_child_path = dir(['C:\Users\noambox\Documents\Sync\Neural data\',parent_path{1}]);
        for child_path = {list_child_path.name}
            if(strfind(child_path{1}, 'L'))
                pathName = ['C:\Users\noambox\Documents\Sync\Neural data\',parent_path{1},'\',child_path{1}];
                fName_table = dir([pathName,'\*.xlsx']);
                for fName = fliplr({fName_table.name})
%                     fName = {'DataRaw_ORI_10Hz_3REP_sub1-.xlsx'}; % DEBUGG
                    if(isempty(strfind(fName{1}, 'Artif')))
                        display(['******* - running on "',parent_path{1},'\',child_path{1},'\',fName{1},'" - *******']);
                        [ data_stats ] = dataStatistics(['\',fName{1}], pathName , data_stats );
                    end
                end              
            end
        end
    end
end