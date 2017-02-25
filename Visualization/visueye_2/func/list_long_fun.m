function [handles, hObject] = list_long_fun(handles, hObject, plot_mode)
%UNTITLED8 Summary of this function goes here
%   plot_mode == 0 - only signal update

if ~ischar(hObject)
    contents = cellstr(get(hObject,'String'));
    selection = contents{get(hObject,'Value')};
else
    selection = hObject;
end

%     [~, tindx] = ismember(handles.table_focus.cell_name, selection);
%     [celd, default_props] = allResponsesGivenCellId(handles, tindx);
    [celd, param, mSig, mStim, mSpikes, handles] = lumpingCellData(handles, selection);

    axes(handles.axes1);
    set(gca,'NextPlot','replacechildren');
    if ~handles.nographics
        showInference_interactive(celd.time, celd.data.Df, celd.data.C, celd.data.S, param);
    end
%% debug section
% celd = buildCellDataStructure(handles, celd, mSig, mStim);
% celd.data.S(isnan(celd.data.S)) = 0;
% h_cpy = handles;
% ori_feature_extraction( celd, h_cpy )
%%
if(plot_mode)
    celd = buildCellDataStructure(handles, celd, mSig, mStim,mSpikes);
    visualizeResponse( handles, celd, param );
end
end

function visualizeResponse(handles, cell_data, param )
    vNames = cell_data.Properties.VariableNames;
    if(any(ismember(cell_data.Properties.VariableNames,'properties')))
        cell_data.props = {cell_data.properties};
    end
    switch cell_data.props{:}.stim_type
       
        case 'ORI'
            %%
            axes(handles.axes1);
            set(gca,'NextPlot','replacechildren');
%             param.axis = [0 cell_data.time(end) min(cell_data.data.C)-eps 1.02.*max(cell_data.data.C)+eps];
            param.scale_s = 1;
            if(any(ismember(cell_data.Properties.VariableNames,'stimulus')))
                param.stimulus = cell_data.stimulus{:};
            else
                param.stimulus = cell_data.stim{:};
            end
            param.stim_type = cell_data.props{1}.stim_type;
            showInference_interactive(cell_data.time, cell_data.data.Df, cell_data.data.C, cell_data.data.S, param)
            
            [log_table,tcd] = oriSummary(handles, cell_data, param); % tcd is tuning curve data
            
            plotOri(handles, tcd , cell_data);
            
            % log
            anova_details = evalc('log_table');
            log_string = [
                'QI: ',num2str(cell_data.data.QI),...
                char(10),'MI: ',num2str(cell_data.data.MI),...
                char(10),'CSD: ',num2str(cell_data.data.CSD),...
                char(10),'*anova p value comparison*',...
                anova_details(193:end)];
            set( handles.analysis_TXT,'String',log_string );
   
        case 'CRP'
            %% 
            axes(handles.axes1);
            set(gca,'NextPlot','replacechildren');
            param.axis = [0 cell_data.time(end) min(cell_data.data.C)-eps 1.3.*max(cell_data.data.C)+eps];
            param.scale_s = 1;
            if(any(ismember(vNames ,'stim')))
                param.stimulus = cell_data.stim{:};
            else
                stim_struct = handles.stim_map(cell2mat(cell_data.fkey));
                param.stimulus = stim_struct.stim;
                lin_s2r = cell2mat(arrayfun(@(a) (handles.samples2remove(1,a):handles.samples2remove(2,a)),(1:size(handles.samples2remove,2)),'UniformOutput',false));
                cell_data.stim{:}.stim(lin_s2r) = [];
                cell_data.stim{:}.partition(:,logical(prod(ismember(cell_data.stim{:}.partition,handles.samples2remove)))) = [];
            end
            param.stim_type = cell_data.props{1}.stim_type;
%             showInference_interactive(cell_data.time, cell_data.data.Df, cell_data.data.C, cell_data.data.S, param)
            
            % summary function
            visual_flag = true;
            [fv_s,fv_c, visu_vec_s ] = crpSummary(handles,cell_data, visual_flag );
            if(~any(fv_s)); return; end;
            
            Table = plotCrp(handles, fv_s, fv_c, visu_vec_s, cell_data); % plots
            
            % log
            tabel_total = Table(:,8);
            tabel_total = cell2mat(tabel_total(3:end)');
            log_string = [
                'On: ',num2str(tabel_total(1)),...
                char(10),'Off: ',num2str(tabel_total(2)),...
                char(10),'Gray: ',num2str(tabel_total(3)),...
                char(10),'Chirp: ',num2str(tabel_total(4)),...
                char(10),'Amp: ',num2str(tabel_total(5))];
            set( handles.analysis_TXT,'String',log_string );
        
        case 'BDN'
            %% 
            bdnSummary(handles,cell_data);
            
        case 'CNT'
            group_label = 0;foo = 0;
            [ mean_response, keys_cpd, max_contrast,cell_data ] = cntSummary(handles,cell_data);   
            plotCnt(handles,mean_response,keys_cpd,max_contrast );
            
            if(handles.plotmore_CHBX.Value)
               sandbox_fn_foldCnt(handles, cell_data);
            end
            
        case 'WGN'
            [kernel,time, all_signals,stimulus] = wgnSummary(handles,cell_data);
            plotWgn( handles,time,kernel,all_signals,stimulus,cell_data )
        
            if(handles.plotmore_CHBX.Value)
%                sandbox_fn_multicell_decode(handles, cell_data);
               sandbox_fn_allcell_decode(handles, cell_data);
%                executeScript(); 
            end
        otherwise
            error('stimulus type is not valid')
    end
    end
