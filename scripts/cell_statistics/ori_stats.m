% for thresh_1 data
% script to find how how many OS DS ON OFF cells
load('C:\Users\noambox\Documents\Sync\Neural data\stim_features\ori_feature_mat_thresh1.mat');

% filter out 14 and 15/3
vault_exper = cellfun(@(faulty) cellfun(@isempty, strfind(table2save.Properties.RowNames,faulty)),{'14-Mar','15-Mar','13-Mar','16-Mar'},...
    'UniformOutput',false);
valid_ind = prod(cell2mat(vault_exper),2);
isnanCell = @(c) cell2mat(cellfun(@(x) sum(isnan(x)),c,'UniformOutput',false));
tbl = table2save(logical(valid_ind),:);

%% print stats ON

sec_indx = tbl.on_off == 1 & isnan(tbl.onoff);
a = nnz(sec_indx);disp(['ON Cells ',num2str(a),' out of ', num2str(height(tbl)),'. OVR -(%) ',num2str(100*a/height(tbl))])
b = nnz(sec_indx & ~isnanCell(tbl.os) & isnanCell(tbl.ds));disp(['ON-OS Cells ',num2str(b),' out of ', num2str(nnz(sec_indx)),'. OVR -(%) ',num2str(100*b/nnz(sec_indx))])
c = nnz(sec_indx & ~isnanCell(tbl.ds));disp(['ON-DS Cells ',num2str(c),' out of ', num2str(nnz(sec_indx)),'. OVR -(%) ',num2str(100*c/nnz(sec_indx))])
%% print stats OFF

sec_indx = tbl.on_off == 0 & isnan(tbl.onoff);
a = nnz(sec_indx);disp(['OFF Cells ',num2str(a),' out of ', num2str(height(tbl)),'. OVR -(%) ',num2str(100*a/height(tbl))])
b = nnz(sec_indx & ~isnanCell(tbl.os) & isnanCell(tbl.ds));disp(['OFF-OS Cells ',num2str(b),' out of ', num2str(nnz(sec_indx)),'. OVR -(%) ',num2str(100*b/nnz(sec_indx))])
c = nnz(sec_indx & ~isnanCell(tbl.ds));disp(['OFF-DS Cells ',num2str(c),' out of ', num2str(nnz(sec_indx)),'. OVR -(%) ',num2str(100*c/nnz(sec_indx))])
%% print stats ON-OFF

sec_indx = ~isnan(tbl.onoff);
a = nnz(sec_indx);disp(['ONOFF Cells ',num2str(a),' out of ', num2str(height(tbl)),'. OVR -(%) ',num2str(100*a/height(tbl))])
b = nnz(sec_indx & ~isnanCell(tbl.os) & isnanCell(tbl.ds));disp(['ONOFF-OS Cells ',num2str(b),' out of ', num2str(nnz(sec_indx)),'. OVR -(%) ',num2str(100*b/nnz(sec_indx))])
c = nnz(sec_indx & ~isnanCell(tbl.ds));disp(['ONOFF-DS Cells ',num2str(c),' out of ', num2str(nnz(sec_indx)),'. OVR -(%) ',num2str(100*c/nnz(sec_indx))])
