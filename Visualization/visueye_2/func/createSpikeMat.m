function out = createSpikeMat(struct_data)
% createSMat(struct_data) handles cases where signal lenght is not the same for all
    S_set = {struct_data(:).S};
    max_df_size = max(cellfun(@length, S_set));
    xq = linspace(0,max_df_size-1,max_df_size);
    eqSize = @(s) z1(reshape(interp1(linspace(0,length(s)-1,length(s)), s, xq,'spline'),[],1));
    out = cellfun(eqSize, S_set, 'UniformOutput',false)';
end