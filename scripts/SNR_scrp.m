% find SNR
% load first data table

for k = 1:numel(handles.table.data)
    df = handles.table.data(k).Df;
    s = handles.table.data(k).S;
    peak_f = max(df);
    
%     baseline samples
    fs = 10;
    bins = 5 * fs; % sec * sr
    num_bin = floor( length(df)/ bins );
    df_mat_binned = reshape(df(1:num_bin*bins,:),bins, num_bin, []);
    perc = .20; % the percent that the hist values under are probabely bsl
    std_vec = std(df_mat_binned, 1);
    sort_std_vec = sort(std_vec);
    SNR(k) = peak_f/mean(sort_std_vec(1:5));
end