% parse data script 
% moves across all files from a .xlmx type in the directory folder and subfolder and execute the parse
% data function
global logvar;
logvar = 'List of files with no artifact:';

exclude_list = ['C:\Users\noambox\Documents\Sync\Neural data\160314 - gcamp6s - old\L1 ',...
    'C:\Users\noambox\Documents\Sync\Neural data\160314 - gcamp6s - old\L2',...
    'C:\Users\noambox\Documents\Sync\Neural data\160314 - gcamp6s - old\L1'
    ];

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
                        display([char(10) '******* - running on "',pathName,'\',fName{1},'" - *******' char(10)]);
%                         [ data_stats ] = dataStatistics(['\',fName{1}], pathName , data_stats );
                        if(isempty(strfind(pathName,'160315')) && isempty(strfind(pathName,exclude_list))) % remove and fix
                            parseData_fn( ['\',fName{1}], pathName );
                                disp('foo');
                        end
                    end
                end              
            end
        end
    end
end