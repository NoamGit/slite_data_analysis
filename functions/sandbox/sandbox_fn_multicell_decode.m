function [  ] = sandbox_fn_multicell_decode( h, cell_data )
% for multicellular rate

% combine traces
on_cell = {'Cell_13','Cell_15','Cell_43', 'Cell_46'}; % 280316
off_cell = {'Cell_4','Cell_9','Cell_42'}; % 280316
on_cell = {'Cell_4','Cell_14','Cell_44','Cell_54'}; %020516

cell_struct = on_cell;
on_ind = ismember(h.table_focus.cell_name,cell_struct);
df = h.table_focus(on_ind,:);
temp = df.data;
temp = (reshape([temp.S],[],numel(cell_struct)));
% temp = bsxfun(@rdivide,temp,max(temp,[],1)); %normalize to max = 1
S_all = sum(temp,2);
cell_data.data.S = reshape(S_all,1,[]);

[kernel,time, all_signals,stimulus] = wgnSummary(h,cell_data);
plotWgn( h,time,kernel,all_signals,stimulus,cell_data )

end

