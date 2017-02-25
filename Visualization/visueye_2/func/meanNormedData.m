function [out, multi_signal,multi_stim, multi_S] = meanNormedData(cell_data,scale)
    if(any(ismember(cell_data.Properties.VariableNames,'properties')))
        l = length(cell_data.data.C);
    else
        l = length(cell_data.data.C);
    end
    [all_data_S, all_data_C, all_data_F] = deal(zeros(l,1));
    all_baseline = 0;
    num_reps = numel(cell_data(:,1));
    full_data = cell_data.data(:);
    baseline_data = [full_data(:).mcmc_samples];
    minmax = @(x) [min(x); max(x)];
    C = reshape([full_data(:).C],[],num_reps);
    S = reshape([full_data(:).S],[],num_reps);
    F = reshape([full_data(:).Df],[],num_reps);
    [minmax_c] = minmax(C);
    [~,ii_dominant] = min(minmax_c(1,:));
    [minmax_s] = minmax(S);
    [~,ii_dominant_s] = min(minmax_s(1,:));
    [minmax_f] = minmax(F);
    baseline = reshape([baseline_data(:).Cb],[],num_reps);
    num_inn_reps = size(cell_data.stim{1}.partition,2);
    multi_struct  = cell(num_reps * num_inn_reps,1);
    stimulus = cell_data.stim{:};
    multi_ind = (1:num_inn_reps:numel(multi_struct));
    min_length = min(cellfun(@(c) c(2) - c(1)+1,num2cell(stimulus.partition,1)));
    for k = 1:num_reps
        multi_struct(multi_ind(k):multi_ind(k)+num_inn_reps-1) = ...
                cellfun(@(c) C((c(1):c(1)+min_length-1),k),num2cell(stimulus.partition,1),'UniformOutput',false)';
        multi_struct_S(multi_ind(k):multi_ind(k)+num_inn_reps-1) = ...
                cellfun(@(c) S((c(1):c(1)+min_length-1),k),num2cell(stimulus.partition,1),'UniformOutput',false)';
        all_data_C = all_data_C + z1(C(:,k));
        all_data_S = all_data_S + z1(S(:,k));
        all_data_F = all_data_F + z1(F(:,k));
        all_baseline = all_baseline + mean(baseline(:,k));
    end
    multi_stim = cellfun(@(c) stimulus.stim((c(1):c(1)+min_length-1)),num2cell(stimulus.partition,1),'UniformOutput',false)';
    multi_S = cellfun(@(c) mapminmax_nc(c,minmax_s(1,ii_dominant_s),minmax_s(2,ii_dominant_s)),multi_struct_S,'UniformOutput',false)';
    multi_signal = cellfun(@(c) mapminmax_nc(c,minmax_c(1,ii_dominant),minmax_c(2,ii_dominant)),multi_struct,'UniformOutput',false);

%     multi_signal = multi_struct;

    % back to dominant scale
    out.S = mapminmax_nc(all_data_S'./num_reps,minmax_s(1,ii_dominant),minmax_s(2,ii_dominant));
    out.C = mapminmax_nc(all_data_C'./num_reps,minmax_c(1,ii_dominant),minmax_c(2,ii_dominant));
    out.F = mapminmax_nc(all_data_F'./num_reps,minmax_f(1,ii_dominant),minmax_f(2,ii_dominant));
    out.baseline = mapminmax_nc(all_baseline'./num_reps,minmax_f(1,ii_dominant),minmax_f(2,ii_dominant));
    out.baseline = mean(baseline(:,ii_dominant));
    
    % if scale is specified then scale according to ca signal
    if( nargin >1 )
        [out.C, norm_fun] = mapminmax_nc(out.C,scale(1),scale(2));
        out.S = norm_fun(out.S);
        out.F = norm_fun(out.F);
        out.baseline = norm_fun(out.baseline);
        multi_signal = cellfun(@(c) norm_fun(c),multi_signal,'UniformOutput',false);
    end
    
    out.S = reshape(out.S,[],1);    
    out.C = reshape(out.C,[],1);
    out.F = reshape(out.F,[],1); 
end