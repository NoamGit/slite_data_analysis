function txt = crpInteractive(empt,event_obj)
% Customizes text of data tips
global crp_features_filt;
global lbl;
global crp_table;
global crp_features;

pos = get(event_obj,'Position');
[~, indx] = ismember( [pos(1),pos(2),pos(3) ],[ crp_features_filt.QI,crp_features_filt.MI,crp_features_filt.CSD ],'rows');
location_txt = ['QI: ',num2str(pos(1)),'  MI: ',num2str(pos(2)),'  CSD: ',num2str(pos(3))];
txt = {['File: ',cell2mat( lbl(indx,1) )],['Location: ',cell2mat( lbl(indx,2))],location_txt};