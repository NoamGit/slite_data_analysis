function [ celd ] = resampleCellData(head, all_d,celd)
%resampleCellData(head, all_d) gets a head index and all duration vectors
%and normalizes all tail vectors (not head) according to head and do
%histcounts to check which samples should be removed

    l = numel(all_d);
    tail = ones(l,1);
    tail(head) = 0; % logical index
    tail_d = all_d;
    forloop_range = setdiff((1:l),head);
    for tail_ind = forloop_range
        
        % create time partition vector of tail according to tail_vec
        tail_vec = ([0;cell2mat(tail_d(tail_ind))]);
        tail_vec_cumsum = cumsum(tail_vec);
        tail_nSamples = celd.properties(tail_ind).signalLenght;
        tail_dt = 1/celd.properties(tail_ind).imaging_rate;
        tail_lintime = linspace(0,(tail_nSamples-1) * tail_dt ,tail_nSamples);
        counts_tail = histcounts(tail_vec_cumsum(2:end) ,[tail_lintime, inf])'; % stim to response
        tail_time = cumsum(cellfun(@sum,(mat2cell(tail_vec(2:end),counts_tail,1))));
        if(numel(tail_time) == celd.properties(tail_ind).signalLenght )
            tail_time = [0; tail_time(1:end-1)];
        end
        
        % N = number of sampels differing the head from tail
        N = celd.properties(tail_ind).signalLenght - celd.properties(head).signalLenght;
        
        % for N largest dt's
        [~, I] = sort(diff(tail_time));
        if( all(diff(I((2:N)) == 1)) ) % checks if the smearing of time is uniform
            binsize = floor(celd.properties(tail_ind).signalLenght/N);
            index2mean = randi(binsize,N,1);
            index2mean = index2mean + (0:N-1)' * binsize;
        else
            index2mean = [I((1:N));I((1:N))+1];
        end
        
        % update data
        celd.data(tail_ind).Df = meanAndThrow( celd.data(tail_ind).Df, index2mean );
        celd.data(tail_ind).S = meanAndThrow( celd.data(tail_ind).S, index2mean );
        celd.data(tail_ind).C = meanAndThrow( celd.data(tail_ind).C, index2mean );
    end
end
