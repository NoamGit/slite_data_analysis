function [cell_data] = removeData(cell_data, s2r)
    lin_s2r = cell2mat(arrayfun(@(a) (s2r(1,a):s2r(2,a)),(1:size(s2r,2)),'UniformOutput',false));
    vNames = cell_data.Properties.VariableNames;
    if(any(ismember(vNames,'properties'))); cell_data.properties.signalLenght = cell_data.properties.signalLenght -numel(lin_s2r); end
    if(any(ismember(vNames,'props'))); cell_data.props{:}.signalLenght = cell_data.props{:}.signalLenght -numel(lin_s2r); end
    if(any(ismember(vNames,'time')))
       temp = cell_data.time;
       cell_data.time = [];
       temp(lin_s2r) = [];
       cell_data.time = temp;  
    end
    cell_data.stim{:}.stim(lin_s2r) = [];
    cell_data.stim{:}.partition(:,logical(prod(ismember(cell_data.stim{:}.partition,s2r)))) = [];
    cell_data.data.raw(lin_s2r) = [];
    cell_data.data.artifact(lin_s2r) = [];
    cell_data.data.Df(lin_s2r) = [];
    cell_data.data.C(lin_s2r) = [];
    cell_data.data.S(lin_s2r) = [];
end
