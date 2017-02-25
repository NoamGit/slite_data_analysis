function out = createDfMat(struct_data)
% createDfMat(struct_data) handles cases where signal lenght is not the same for all
    C_set = {struct_data(:).C};
    max_df_size = max(cellfun(@length, C_set));
    xq = linspace(0,max_df_size-1,max_df_size);
    eqSize = @(s) z1(reshape(interp1(linspace(0,length(s)-1,length(s)), s, xq,'spline'),[],1));
    out = cellfun(eqSize, C_set, 'UniformOutput',false)';
end

