function [ kernel, time ,sta_container, stimulus ] = wgnSummary( handles,cell_data, param  )
%wgnSummary( handles,cell_data, param  )returns a 1D vector time kernel
%according to reverse correlation (STA) technique
%   1D case of BDN
% constants
rf_bounds = 0; % msec
stim_fs = 30; % Hz
num_rep = 2;
off_set = 200;

% define parameters
if(str2double(get(handles.param_ETXT,'String')))
    peri_stim_time = str2double(get(handles.param_ETXT,'String'))+rf_bounds; % in msec
else
    peri_stim_time = 500; % in msec
end
off_set_samples = floor(off_set / (1e3/cell_data.props{1}.imaging_rate));
pstl = floor(peri_stim_time / (1e3/cell_data.props{1}.imaging_rate)); % peri-stimulus latency in samples
l = cell_data.props{1}.signalLenght;
dt = 1/cell_data.props{1}.imaging_rate;
if(cell_data.props{:}.isnew)
    stimulus = handles.wgn_src.new;
else
    stimulus = handles.wgn_src.old;
end

% find all relevant samples  
resampling_ratio = cell_data.props{:}.imaging_rate/stim_fs;
w = cell_data.data.S(cell_data.data.S ~= 0);
w = w/max(w); % norma weights
spike_indx_re = find(cell_data.data.S ~= 0);

% map them to stimulus
spike_indx_stim = spike_indx_re./resampling_ratio;
stimulus = repmat(stimulus,2,1); 
stimulus = stimulus(:);
tidx_2 = repmat(spike_indx_stim,(pstl+off_set_samples)/resampling_ratio + 1,1);
tidx_3 = [fliplr((0:off_set_samples/resampling_ratio))'; -(1:pstl/resampling_ratio)']; %concatenating offset and pstl
tidx = bsxfun(@plus, tidx_2, tidx_3);
tidx_remove = find(tidx(1,:) > length(stimulus) | tidx(end,:) < 1); % remove indx out of range
tidx(:, tidx_remove) = [];
w(:, tidx_remove) = [];

% resampling_ratio = cell_data.props{:}.imaging_rate/stim_fs;
% resampling_ratio
% rs_vec = repmat((1:ceil(length(stimulus)/resampling_ratio)),2,1);
% rs_vec = repmat((1:ceil(length(stimulus)/resampling_ratio)),2,1);
% rs_vec = rs_vec(:);
% rs_vec = rs_vec(1:length(stimulus));

% build the 1D container
sta_container = zeros(size(tidx));
for k = 1:size(tidx,2)
    try
        sta_container(:,k) = w(k) .* stimulus(tidx(:,k));
    catch exception
        disp((tidx(:,k)));
    end
end
if(isempty(sta_container))
    kernel = zeros(1,size(tidx,1));
else
    kernel = mean(sta_container,2);
end
sta_container = bsxfun(@minus,sta_container ,kernel);
kernel = kernel - mean(kernel);
time = tidx_3 / stim_fs;

end

