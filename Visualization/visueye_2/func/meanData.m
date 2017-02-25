function [out, multi_signal] = meanData(cell_data)
    if(any(ismember(cell_data.Properties.VariableNames,'properties')))
        l = length(cell_data.data.C);
    else
        l = length(cell_data.data.C);
    end
    [all_data_S, all_data_C, all_data_F] = deal(zeros(1,l));
    all_baseline = 0;
    num_reps = numel(cell_data(:,1));
    num_inn_reps = size(cell_data.stim{1}.partition,2);
    multi_struct  = cell(num_reps * num_inn_reps,1);
    multi_ind = (1:num_inn_reps:numel(multi_struct));
    for k = 1:num_reps
        all_data_C = all_data_C + cell_data(k,:).data.C;
        all_data_S = all_data_S + cell_data(k,:).data.S;
        all_data_F = all_data_F + cell_data(k,:).data.Df';
        all_baseline = all_baseline + mean(cell_data(k,:).data.mcmc_samples.Cb);
    end
    out.S = reshape(all_data_S./num_reps,[],1);
    out.C = reshape(all_data_C./num_reps,[],1);
    out.F = reshape(all_data_F./num_reps,[],1);
    out.baseline = all_baseline./num_reps;
end
