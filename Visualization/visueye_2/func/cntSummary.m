function [ mean_response, keys_cpd, max_contrast, cell_data, group ] = cntSummary( handles,cell_data, param )
%cntSummary( handles,cell_data, param ) 
%   accumulates responses according to number of bins

% leave it: find cpd's for old stimulus...
% TODO: 1. fold response
%       2. fix plots according to WebVision

% quantize each row according to experiment value
if ~ismember('stimulus', cell_data.Properties.VariableNames)
    cell_data(:,'stimulus') = {cell_data.stim};
    cell_data(:,'properties') = cell_data.props;
end
contrast_vec = cell_data.stimulus{:}.stim(:,1);
contrast_vec = quantiz(contrast_vec,linspace(1,max(contrast_vec),6))+1;
map_vec = zeros(max(contrast_vec),1); map_vec(unique(contrast_vec)) = (1:numel(unique(contrast_vec)));
contrast_vec = map_vec(contrast_vec);
cpd_vec = cell_data.stimulus{:}.stim(:,2);
cpd_vec = quantiz(cpd_vec,[7 8 11 14 26 58 70],(1:7))+1;
map_vec = zeros(max(cpd_vec),1); map_vec(unique(cpd_vec)) = (1:numel(unique(cpd_vec)));
cpd_vec = map_vec(cpd_vec);

if(cell_data.properties.isnew)
    keys_cpd = [0.04, 0.06, 0.08 0.1 0.2 0.4 0.6]; % cpds
    cpd_dict = containers.Map(keys_cpd,unique(cpd_vec));
    keys_dict = [0 0.1 0.3 0.5 0.7 0.9];
    contrast_dict = containers.Map(keys_dict,unique(contrast_vec));
else
    error('cntSummary:: FIND OLD CPD''s');
end

% accumulate for each 2D index
group = [contrast_vec cpd_vec+max(contrast_vec)];

% % stimulus includes off frame with on frame
% on_time = 2; % sec
% subs_shifted_1 = circshift(subs(:,1),on_time * cell_data.properties.imaging_rate);
% % subs_shifted_2  = medfilt1(foo,3 * on_time * cell_data.properties.imaging_rate);
% foo = subs(:,1);
% foo( subs(:,1) == 1 ) = subs_shifted_1( subs(:,1) == 1 );
% subs_shifted_1 = medfilt1(foo,3);
% subs(:,1) = subs_shifted_1;
% cell_data.stim{:}.stim(:,1) = subs(:,1).^2;

mean_response_full = accumarray(group, cell_data.data.S,[],@mean);
mean_response = mean_response_full(:,(min(group(:,2)):end));
% disp(mean_response);
[~,ii] = max(mean_response,[],1);
max_contrast = keys_dict(ii);
end

