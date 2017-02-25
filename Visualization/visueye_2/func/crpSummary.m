function [ feature_vec_S, feature_vec_C , visu_vec_s] = crpSummary( handles,cell_data, visual_flag )
%crpSummary( handles,cell_data, param ) yield a 20 dim representation of
%the features vector [on off gray chirp(1,5) amp(1,12)]

    % constants
    crp_part = [13    20    26    60]-2; % 5 states extracted with crp_partition_scrp
    amp_part = [16    27    35    48    57    71    79    90   101   112   122]-2; % 11 states extracted with crp_partition_scrp

    stim_struct = handles.stim_map(cell2mat(cell_data.fkey));
    lin_s2r = cell2mat(arrayfun(@(a) (handles.samples2remove(1,a):handles.samples2remove(2,a)),(1:size(handles.samples2remove,2)),'UniformOutput',false));
    stim_struct.stim(lin_s2r) = [];
    stim_struct.partition(:,logical(prod(ismember(stim_struct.partition,handles.samples2remove)))) = [];
    i_crp = find(stim_struct.stim > 1.5 & stim_struct.stim < 2.5);
    i_amp = find(stim_struct.stim > 3 & stim_struct.stim < 4);
    i_on = find(stim_struct.stim == 1);
    i_off = find(stim_struct.stim == 0);
    i_gray = find(stim_struct.stim > 0.49 & stim_struct.stim < 0.51);
    
    % fold signal and extract num of spikes from every state
    on_C_mu = mean(cell_data.data.C(i_on));
    on_S_mu = mean(cell_data.data.S(i_on));
    off_C_mu = mean(cell_data.data.C(i_off));
    off_S_mu = mean(cell_data.data.S(i_off));
    gray_C_mu =  mean(cell_data.data.C(i_gray));
    gray_S_mu = mean(cell_data.data.S(i_gray));
    
    on_C_std = std(cell_data.data.C(i_on));
    on_S_std = std(cell_data.data.S(i_on));
    off_C_std = std(cell_data.data.C(i_off));
    off_S_std = std(cell_data.data.S(i_off));
    gray_C_std =  std(cell_data.data.C(i_gray));
    gray_S_std = std(cell_data.data.S(i_gray));
    
    const_vec = cell(length(i_on) + length(i_off)+length(i_gray),2);
    const_vec(:,1) = num2cell([cell_data.data.S(i_on) cell_data.data.S(i_off) cell_data.data.S(i_gray)]);
    const_vec(:,2) = cat(1,num2cell(repmat('on',length(i_on),1),2) ,num2cell(repmat('off',length(i_off),1),2),num2cell(repmat('gray',length(i_gray),1),2));
    
    [crp_S, crp_C, amp_S, amp_C] = deal(cell(size(stim_struct.partition,2),1));
    for n = 1:size(stim_struct.partition,2)
        i_part1 = find(i_crp >= stim_struct.partition(1,n) & i_crp <= stim_struct.partition(2,n));
        crp_part_2 = [crp_part numel(i_part1)];
        crp_part_1 = [1 crp_part_2(1:end-1)+1];
        crp_C{n} = arrayfun(@(x,y) cell_data.data.C(i_part1(x:y)),crp_part_1,crp_part_2,'UniformOutput',false);
        crp_S{n} = arrayfun(@(x,y) cell_data.data.S(i_part1(x:y)),crp_part_1,crp_part_2,'UniformOutput',false);

        i_part2 = find(i_amp >= stim_struct.partition(1,n) & i_amp <= stim_struct.partition(2,n));
        amp_part_2 = [amp_part numel(i_part2)];
        amp_part_1 = [1 amp_part_2(1:end-1)+1];
        amp_C{n} = arrayfun(@(x,y) cell_data.data.C(i_part2(x:y)),amp_part_1,amp_part_2,'UniformOutput',false);
        amp_S{n} = arrayfun(@(x,y) cell_data.data.S(i_part2(x:y)),amp_part_1,amp_part_2,'UniformOutput',false);
    end
    
    crp_S_mat = cell2mat([crp_S{:}])';
    amp_S_mat = cell2mat([amp_S{:}])';
    visu_vec_s = cell(size(const_vec,1) + numel(crp_S_mat) + numel(amp_S_mat),2);
    visu_vec_s(:,1) = [const_vec(:,1);num2cell(crp_S_mat);num2cell(amp_S_mat)];
    visu_vec_s(:,2) = cat(1,const_vec(:,2),num2cell(repmat('chirp',length(crp_S_mat),1),2),...
                        num2cell(repmat('amp',length(amp_S_mat),1),2));
                    
    crp_C_mu = mean(cell2mat(cellfun(@(x) cellfun(@mean,x), crp_C,'UniformOutput',false)),1);
    crp_S_mu = mean(cell2mat(cellfun(@(x) cellfun(@mean,x), crp_S,'UniformOutput',false)),1);
    amp_C_mu = mean(cell2mat(cellfun(@(x) cellfun(@mean,x), amp_C,'UniformOutput',false)),1);
    amp_S_mu = mean(cell2mat(cellfun(@(x) cellfun(@mean,x), amp_S,'UniformOutput',false)),1);
    
    crp_C_std = mean(cell2mat(cellfun(@(x) cellfun(@std,x), crp_C,'UniformOutput',false)),1);
    crp_S_std = mean(cell2mat(cellfun(@(x) cellfun(@std,x), crp_S,'UniformOutput',false)),1);
    amp_C_std = mean(cell2mat(cellfun(@(x) cellfun(@std,x), amp_C,'UniformOutput',false)),1);
    amp_S_std = mean(cell2mat(cellfun(@(x) cellfun(@std,x), amp_S,'UniformOutput',false)),1);
    
    % return feature vector
    feature_mu_S = z1([on_S_mu off_S_mu gray_S_mu crp_S_mu amp_S_mu]);
    feature_vec_S = [feature_mu_S;[ on_S_std  off_S_std  gray_S_std  crp_S_std  amp_S_std]];
    feature_vec_S = feature_vec_S(:);
    
    feature_mu_C = z1([on_C_mu off_C_mu gray_C_mu crp_C_mu amp_C_mu]);
    feature_vec_C = [feature_mu_C;[ on_C_std  off_C_std  gray_C_std  crp_C_std  amp_C_std]];
    feature_vec_C = feature_vec_C(:);
    if visual_flag
        feature_vec_S = feature_mu_S;
        feature_vec_C = feature_mu_C;
    end
end

