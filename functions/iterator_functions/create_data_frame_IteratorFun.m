function [ tbl_out ] = create_data_frame_IteratorFun( celd, h )
% function to create the data structre for the revised 'chirp' analysis discussed
% with Shy on the 16/2/17

global celltype
global stim_map
alpha = 0.05; % confidence value for OnOff test
size_on_off_wn = 2; % 0.2 sec win

% filter out vaulty experiments
vault_exper = cellfun(@(faulty) (strfind(celd.cell_id,faulty)),{'14-Mar','15-Mar','13-Mar','16-Mar'});
if any(~cellfun(@isempty,vault_exper))
    tbl_out = [];
    return
end

%  finding gcamp type
if(any( ismember(celd.Properties.VariableNames,'properties') ))
    gcamp_type = {celd.properties.gcampType};
    key = cell2mat(celd.fkey);
    stim = stim_map(key).stim;
    stim_partition = stim_map(key).partition;
elseif(any( ismember(celd.Properties.VariableNames,'props') ))
    gcamp_type = {celd.props{:}.gcampType};
    stim = celd.stim{:}.stim;
    stim_partition = celd.stim{:}.partition  ;
end

% find cell_id and update type
[~,ind] = ismember(celd.cell_id, celltype.Properties.RowNames);
cell_type = {' '};
if(ind)
    switch celltype(ind,:).ONorOFF
        case 1
            cell_type = {'on'};
        case 2
            cell_type = {'off'};
   end
    if ~isnan(celltype(ind,:).ONorOFF)
        cell_type = {'on_off'};
    end
end

% calculate mean response of each state for crp
i_crp = find(stim > 1.5 & stim < 2.5);
i_amp = find(stim > 3 & stim < 4);
i_gray = find(stim > 0.49 & stim < 0.51);
i_off = find(stim < 0.01);
i_on = find(stim > 0.51 & stim < 1.01);
mat_s = cell2mat(celd.multi_data.multi_spikes');
mat_c = cell2mat(celd.multi_data.multi_signal');

%  creating the chirp and amp matrix
crp_ind = find(i_crp >= stim_partition(1,1) & i_crp <= stim_partition(2,1));
crp_mat_s = {mat_s(i_crp(i_crp(crp_ind)< length(mat_s)),:)'};
crp_mat_c = {mat_c(i_crp(i_crp(crp_ind)< length(mat_s)),:)'};

amp_ind = find(i_amp >= stim_partition(1,1) & i_amp <= stim_partition(2,1));
amp_mat_s = {mat_s(i_amp(i_amp(amp_ind)< length(mat_s)),:)'};
amp_mat_c = {mat_c(i_amp(i_amp(amp_ind)< length(mat_s)),:)'};

%  lumping all other parts
on_ind = i_on >= stim_partition(1,1) & i_on  <= stim_partition(2,1);
off_ind = i_off >= stim_partition(1,1) & i_off <= stim_partition(2,1);
gray_ind = i_gray >= stim_partition(1,1) & i_gray  <= stim_partition(2,1);
on_mat_s = {mat_s(i_on(i_on(on_ind)< length(mat_s)),:)'};
off_mat_s = {mat_s(i_off(i_off(off_ind) < length(mat_s)),:)'};
gray_mat_s = {mat_s(i_gray(i_gray(gray_ind)< length(mat_s)),:)'};

%  do elementary ON or OFF type cell test in case there is no cell_type
if(strcmp(cell_type,' '))
    
    % on or off test
    [p_on_off,~, stats] = anova1([on_mat_s{:}(:); off_mat_s{:}(:)],...
        [repmat({'on'},size(on_mat_s{:}(:)));  repmat({'off'},size(off_mat_s{:}(:))) ],'off');
%     figure(6);foo = mat_s; plot([mat_s' ;celd.stim{:}.stim(1:length(foo))']');
    if p_on_off < alpha
        [~,argmax] = max(stats.means);
        cell_type = stats.gnames(argmax);
        disp([char(10),' #new ',cell_type{:},' cell!',char(10)]);
    end
    
    % on-off test
    on_off_values = [ reshape(on_mat_s{:}(:,1:size_on_off_wn),[],1); reshape(off_mat_s{:}(:,1:size_on_off_wn),[],1)];
    tail_values = [ reshape(on_mat_s{:}(:,size_on_off_wn+1:end),[],1); reshape(off_mat_s{:}(:,size_on_off_wn+1:end),[],1)];
    [p_onoff_1,~, stats_on_is_off] = anova1(on_off_values,[repmat({'on'},numel(on_mat_s{:}(:,1:size_on_off_wn) ),1);...
        repmat({'off'},numel(off_mat_s{:}(:,1:size_on_off_wn)),1) ],'off');
    [p_onoff_2,~, stats_no_on_off] = anova1([ on_off_values; tail_values ],...
        [ repmat({'short'},length(on_off_values),1); repmat({'tail'},length(tail_values),1)],'off');
    if p_onoff_1 >= alpha && p_onoff_2 < alpha && mean(stats_on_is_off.means) > mean(stats_no_on_off.means)
        cell_type = {'on_off'};
        disp([char(10),' @ new on-off!',char(10)]);
    end
end

%  combine into tble_output
cell_id = celd.cell_id;
tbl_out = table(cell_type, gcamp_type,on_mat_s,off_mat_s,gray_mat_s,crp_mat_s ,crp_mat_c,amp_mat_s,amp_mat_c,'RowNames',cell_id);

end

