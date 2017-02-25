function stim_out = quantizeOriStim(cell_data)
    if isstruct(cell_data.stim{:})
        stim_in = cell_data.stim{:}.stim;
    else
        stim_in = cell_data.stim{:};
    end
    dt = 1/cell_data.props{1}.imaging_rate;
    time = linspace(0,(cell_data.props{1}.signalLenght-1) * dt,cell_data.props{1}.signalLenght);
    num_val_stim = 9;
    [partition, ~] = lloyds(stim_in,num_val_stim);
    [~,stim] = quantiz(stim_in,partition,(1:num_val_stim));
    rise_samp = logical([0 diff(stim) >0]);
    drop_samp = logical([0 diff(stim) < 0]);
    stim(rise_samp) = stim(logical([0 0 0 rise_samp(3:end-3)]));
    stim(drop_samp) = stim(logical([drop_samp(4:end) 0 0 0]));
    stim_out = stim;
end