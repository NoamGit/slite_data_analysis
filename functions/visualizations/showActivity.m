function [  ] = showActivity( filePath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

M = load(filePath);
df = M.dataframe; clear M;

% TODO: insert stimulus times
% stim = loadStimVetor(df.properties.stim_type)

c_mat = varfun(@(x) (x.C), df.cell_table);
df_mat = varfun(@(x) (x.Df), df.cell_table);
matplot_examples = (1:df.properties.numCell);
c_mat_array = reshape(table2array(c_mat),df.properties.signalLenght,[]);
df_mat_array = reshape(table2array(df_mat),df.properties.signalLenght,[]);
figure(9); 
subplot(121);imagesc(c_mat_array(:,matplot_examples)');
title('Calcium estimation matrix');ylabel('time[sec]');
subplot(122);imagesc(df_mat_array(:,matplot_examples)');
title('fluo_{raw} matrix');ylabel('time[sec]')

s_mat_norm = reshape(table2array(varfun(@(x) z1(x.S), df.cell_table)), df.properties.signalLenght,[]);
s_mat_norm = sparse(s_mat_norm);
s_mat = reshape(table2array(varfun(@(x) x.S, df.cell_table)), df.properties.signalLenght,[]);
s_mat = sparse(s_mat);

figure(10);
subplot(122);spyc(s_mat_norm(:,matplot_examples)','Parula',1);
daspect([170 2 1])
title('Raster plot');ylabel('cell');xlabel('time[sec]')
subplot(121);imagesc(c_mat_array(:,matplot_examples)');
colorbar;
daspect([170 2 1])
title('Calcium estimation matrix');ylabel('time [sec]')

end

function [ out ] = loadStimVetor(type)
    switch type
        case 'ORI'
            
        case 'CNT'
            
        case 'CRP'
            
        otherwise
            out = 0;
    end
end